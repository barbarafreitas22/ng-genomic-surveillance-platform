import logging
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import configparser

import pandas as pd
from Bio import SeqIO

logger = logging.getLogger(__name__)

try:
    import psutil as _psutil
except ImportError:
    _psutil = None

def _spades_resources():
    threads = min(os.cpu_count() or 4, 8)
    if _psutil is not None:
        available_gb = _psutil.virtual_memory().available / 1e9
        memory = min(int(available_gb * 0.35), 32)
        return str(threads), str(max(memory, 4))
    return str(threads), "4"

from .paths import PROJECT_ROOT, CONFIG_PATH, REFERENCE_FA1090
from .models import AssemblyResult, AssemblyStats

config = configparser.ConfigParser()
config.read(CONFIG_PATH)
paths = config["PATHS"]

SPADES_PATH = paths.get("SPADES", "spades.py")
QUAST_PATH  = paths.get("QUAST", "quast.py")


def run_cmd(cmd):
    process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    out, err = process.communicate()
    return process.returncode, out, err


def _run_spades(r1, r2, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    _threads, _memory = _spades_resources()
    spades_out = outdir / "spades_output"
    read_args = ["--pe1-1", str(r1), "--pe1-2", str(r2)] if r2 else ["--s1", str(r1)]
    label     = "SPAdes PE" if r2 else "SPAdes SE"
    cmd = [
        SPADES_PATH, *read_args,
        "-o",        str(spades_out),
        "--threads", _threads,
        "--memory",  _memory,
        "--only-assembler",
        "-k", "33,55,77",
        "--cov-cutoff", "auto",
    ]
    rc, out, err = run_cmd(cmd)
    logs = f"[{label}]\n{out}\n{err}\n"
    if rc != 0:
        if rc in (-9, 137):
            logs += f"\n[ERROR] SPAdes killed by OOM (rc={rc}) — cgroup memory limit reached.\n"
        return None, logs
    contigs = spades_out / "contigs.fasta"
    if not contigs.exists():
        logs += "\n[ERROR] contigs.fasta not found.\n"
        return None, logs
    return contigs, logs


def _run_plasmid_spades(r1, r2, outdir, pre_trimmed: bool = False):
    outdir.mkdir(parents=True, exist_ok=True)
    _threads, _memory = _spades_resources()
    spades_out = outdir / "plasmid_spades_output"
    read_args = ["--pe1-1", str(r1), "--pe1-2", str(r2)] if r2 else ["--s1", str(r1)]
    label     = "plasmidSPAdes PE" if r2 else "plasmidSPAdes SE"
    cmd = [
        SPADES_PATH, "--plasmid", *read_args,
        "-o",        str(spades_out),
        "--threads", _threads,
        "--memory",  _memory,
    ]
    if pre_trimmed:
        cmd += ["--only-assembler"]
    rc, out, err = run_cmd(cmd)
    logs = f"[{label}]\n{out}\n{err}\n"
    if rc != 0:
        return None, logs
    contigs = spades_out / "contigs.fasta"
    if not contigs.exists():
        logs += "\n[INFO] plasmidSPAdes produced no output (no plasmids detected).\n"
        return None, logs
    return contigs, logs


def _has_plasmid_content(contigs_path: Path, min_len: int = 500) -> bool:
    """Return True if the plasmid assembly has at least one contig >= min_len bp."""
    if not contigs_path or not contigs_path.exists():
        return False
    return any(len(record.seq) >= min_len for record in SeqIO.parse(str(contigs_path), "fasta"))


def run_assembly_pipeline(sample_id, r1, r2=None, outdir=None, run_plasmid: bool = False,
                          pre_trimmed: bool = False):
    if outdir is None:
        outdir = PROJECT_ROOT / "results" / "assembly" / sample_id
    else:
        outdir = Path(outdir)

    logger.info("Assembly | %s | start pe=%s", sample_id, r2 is not None)

    contigs, logs = _run_spades(r1, r2, outdir)

    if contigs and contigs.exists():
        canonical = outdir / f"{sample_id}.fasta"
        if contigs != canonical:
            shutil.copy2(contigs, canonical)
            tool_dir = contigs.parent
            contigs = canonical
            shutil.rmtree(tool_dir, ignore_errors=True)

    plasmid_contigs = None
    ref = REFERENCE_FA1090 if REFERENCE_FA1090.exists() else None

    with ThreadPoolExecutor(max_workers=2) as _pool:
        quast_fut = (
            _pool.submit(run_quast, contigs, outdir, ref)
            if contigs and contigs.exists()
            else None
        )
        plasmid_fut = None
        if run_plasmid:
            plasmid_fut = _pool.submit(_run_plasmid_spades, r1, r2, outdir, pre_trimmed)

        if quast_fut is not None:
            _, qlogs = quast_fut.result()
            logs += qlogs

        if plasmid_fut is not None:
            _pc, _pl = plasmid_fut.result()
            logs += _pl
            if _has_plasmid_content(_pc):
                plasmid_contigs = _pc
            elif _pc and _pc.parent.exists():
                shutil.rmtree(_pc.parent, ignore_errors=True)

    if plasmid_contigs and plasmid_contigs.exists():
        canonical_plasmid = outdir / f"{sample_id}_plasmid.fasta"
        if plasmid_contigs != canonical_plasmid:
            shutil.copy2(plasmid_contigs, canonical_plasmid)
            shutil.rmtree(plasmid_contigs.parent, ignore_errors=True)
            plasmid_contigs = canonical_plasmid

    if contigs and contigs.exists():
        logger.info("Assembly | %s | done contigs=%s", sample_id, contigs.name)
    else:
        logger.error("Assembly | %s | failed — no contigs produced", sample_id)

    return AssemblyResult(
        contigs=contigs,
        logs=logs,
        plasmid_contigs=plasmid_contigs,
    )


ASSEMBLY_QC_THRESHOLDS = {
    "min_genome_fraction": 90.0,
    "min_n50":             25_000,
    "max_contigs":           300,
    "min_total_len":    1_900_000,
    "max_total_len":    2_500_000,
    "min_gc":               50.0,
    "max_gc":               56.0,
}


def parse_quast_report(quast_dir: Path) -> dict:
    report = quast_dir / "report.tsv"
    if not report.exists():
        return {}
    series = pd.read_csv(report, sep="\t", header=None, index_col=0).iloc[:, 0]
    _int   = {"# misassemblies": "misassemblies", "# misassembled contigs": "misassembled_contigs", "NGA50": "nga50"}
    _float = {"Genome fraction (%)": "genome_fraction", "Duplication ratio": "duplication_ratio",
               "# mismatches per 100 kbp": "mismatches_per_100kbp", "# indels per 100 kbp": "indels_per_100kbp"}
    out: dict = {}
    for src, dst in _int.items():
        if src in series.index:
            try: out[dst] = int(float(series[src]))
            except (ValueError, TypeError): pass
    for src, dst in _float.items():
        if src in series.index:
            try: out[dst] = float(series[src])
            except (ValueError, TypeError): pass
    return out


def evaluate_assembly_qc(stats: dict) -> dict:
    """
    Evaluate assembly stats against ASSEMBLY_QC_THRESHOLDS.

    Returns:
        {
          "status": "pass" | "warn" | "fail",
          "flags":  [list of human-readable failure reasons],
        }

    "warn" = soft thresholds only (N50, contigs); "fail" = hard thresholds
    (core genes, total length, GC).
    """
    t = ASSEMBLY_QC_THRESHOLDS
    flags: list[str] = []
    hard_fail = False

    gf = stats.get("completeness")
    if gf is not None and gf < t["min_genome_fraction"]:
        flags.append(f"Core genes {gf:.1f}% < {t['min_genome_fraction']}%")
        hard_fail = True

    tl = stats.get("total_len")
    if tl is not None:
        if tl < t["min_total_len"]:
            flags.append(f"Assembly too short ({tl/1e6:.2f} Mb < {t['min_total_len']/1e6:.2f} Mb)")
            hard_fail = True
        elif tl > t["max_total_len"]:
            flags.append(f"Assembly too long ({tl/1e6:.2f} Mb > {t['max_total_len']/1e6:.2f} Mb — possible contamination)")
            hard_fail = True

    gc = stats.get("gc_pct")
    if gc is not None and not (t["min_gc"] <= gc <= t["max_gc"]):
        flags.append(f"GC% {gc:.1f}% outside expected range {t['min_gc']}–{t['max_gc']}%")
        hard_fail = True

    n50 = stats.get("n50")
    if n50 is not None and n50 < t["min_n50"]:
        flags.append(f"N50 {n50/1e3:.1f} kb < {t['min_n50']/1e3:.0f} kb (fragmented)")

    nc = stats.get("n_contigs")
    if nc is not None and nc > t["max_contigs"]:
        flags.append(f"{nc} contigs > {t['max_contigs']} (fragmented)")

    if hard_fail:
        status = "fail"
    elif flags:
        status = "warn"
    else:
        status = "pass"

    return {"status": status, "flags": flags}


def parse_contigs_stats(path: Path) -> AssemblyStats | None:
    if not path or not path.exists():
        return None

    lengths, gc_count, total_bases = [], 0, 0
    for record in SeqIO.parse(str(path), "fasta"):
        s = str(record.seq).upper()
        lengths.append(len(s))
        gc_count += s.count("G") + s.count("C")
        total_bases += len(s)

    if not lengths:
        return None

    lengths.sort(reverse=True)
    cum, n50 = 0, 0
    for length in lengths:
        cum += length
        if cum >= total_bases / 2:
            n50 = length
            break

    _stats_d: dict = {
        "n_contigs":   len(lengths),
        "total_len":   total_bases,
        "n50":         n50,
        "largest":     lengths[0],
        "gc_pct":      gc_count / total_bases * 100 if total_bases else 0.0,
        "contigs_500": sum(1 for length in lengths if length >= 500),
        "completeness": 0.0,
    }

    quast_dir = path.parent.parent / "quast"
    quast = parse_quast_report(quast_dir)
    for k in ("misassemblies", "misassembled_contigs", "nga50",
              "duplication_ratio", "mismatches_per_100kbp", "indels_per_100kbp"):
        if k in quast:
            _stats_d[k] = quast[k]

    try:
        from backend.phylogeny.cgmlst import cached_core_genome_completeness
        core = cached_core_genome_completeness(str(path))
    except Exception:
        core = None
    if core:
        _stats_d["completeness"]     = core["core_completeness_pct"]
        _stats_d["core_genes_found"] = core["core_genes_found"]
        _stats_d["core_genes_total"] = core["core_genes_total"]

    qc = evaluate_assembly_qc(_stats_d)
    return AssemblyStats(
        n_contigs=_stats_d["n_contigs"],
        total_len=_stats_d["total_len"],
        n50=_stats_d["n50"],
        largest=_stats_d["largest"],
        gc_pct=_stats_d["gc_pct"],
        contigs_500=_stats_d["contigs_500"],
        completeness=_stats_d["completeness"],
        qc_status=qc["status"],
        qc_flags=qc["flags"],
        core_genes_found=_stats_d.get("core_genes_found"),
        core_genes_total=_stats_d.get("core_genes_total"),
        misassemblies=_stats_d.get("misassemblies"),
        misassembled_contigs=_stats_d.get("misassembled_contigs"),
        nga50=_stats_d.get("nga50"),
        duplication_ratio=_stats_d.get("duplication_ratio"),
        mismatches_per_100kbp=_stats_d.get("mismatches_per_100kbp"),
        indels_per_100kbp=_stats_d.get("indels_per_100kbp"),
    )


def run_quast(contigs, outdir, reference=None):
    quast_out = outdir / "quast"
    quast_out.mkdir(exist_ok=True)

    cmd = [QUAST_PATH, str(contigs), "-o", str(quast_out)]
    if reference and Path(reference).exists():
        cmd += ["-r", str(reference)]

    _, out, err = run_cmd(cmd)
    logs = f"[QUAST]\n{out}\n{err}\n"

    return quast_out, logs
