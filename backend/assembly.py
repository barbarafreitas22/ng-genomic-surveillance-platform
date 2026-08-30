import logging
import os
import shutil
import subprocess
from pathlib import Path
import configparser

from Bio import SeqIO

logger = logging.getLogger(__name__)

try:
    import psutil as _psutil
except ImportError:
    _psutil = None


def _cgroup_memory_gb() -> float | None:
    v2 = Path("/sys/fs/cgroup/memory.max")
    if v2.exists():
        value = v2.read_text().strip()
        if value != "max":
            try:
                return int(value) / (1024 ** 3)
            except ValueError:
                pass

    v1 = Path("/sys/fs/cgroup/memory/memory.limit_in_bytes")
    if v1.exists():
        try:
            value = int(v1.read_text().strip())
            if value < (1 << 60):
                return value / (1024 ** 3)
        except ValueError:
            pass

    return None


def _spades_resources():
    threads = min(os.cpu_count() or 4, 8)
    container_mem = _cgroup_memory_gb()

    if container_mem is not None:
        memory = int(container_mem * 0.70)
        memory = max(memory, 2)
        memory = min(memory, 8)
        return str(threads), str(memory)

    if _psutil is not None:
        available_gb = _psutil.virtual_memory().available / (1024 ** 3)
        memory = int(available_gb * 0.50)
        memory = max(memory, 2)
        memory = min(memory, 8)
        return str(threads), str(memory)

    return str(threads), "4"


from .paths import CONFIG_PATH, RESULTS_DIR
from .models import AssemblyResult, AssemblyStats

config = configparser.ConfigParser()
config.read(CONFIG_PATH)
paths = config["PATHS"]

SPADES_PATH = paths.get("SPADES", "spades.py")


def _run_spades(r1, r2, outdir):
    outdir.mkdir(parents=True, exist_ok=True)

    threads, memory = _spades_resources()

    spades_out = outdir / "spades_output"
    log_path = outdir / "spades.log"

    if r2:
        read_args = [
            "--pe1-1", str(r1),
            "--pe1-2", str(r2),
        ]
        label = "SPAdes PE"
    else:
        read_args = [
            "--s1", str(r1),
        ]
        label = "SPAdes SE"

    cmd = [
        SPADES_PATH,
        "--isolate",
        *read_args,
        "-o", str(spades_out),
        "--threads", threads,
        "--memory", memory,
    ]

    logger.info(
        "Assembly | %s | SPAdes threads=%s memory=%sGB",
        label,
        threads,
        memory,
    )

    with open(log_path, "w", encoding="utf-8") as log:
        process = subprocess.run(
            cmd,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )

    logs = log_path.read_text(encoding="utf-8", errors="replace")

    if process.returncode != 0:
        if process.returncode in (-9, 137):
            logs += (
                "\n[ERROR] SPAdes killed by OOM "
                f"(rc={process.returncode}) — cgroup memory limit reached.\n"
            )
        else:
            logs += (
                f"\n[ERROR] SPAdes failed (rc={process.returncode}).\n"
            )
        return None, f"[{label}]\n{logs}"

    contigs = spades_out / "contigs.fasta"

    if not contigs.exists():
        logs += "\n[ERROR] contigs.fasta not found.\n"
        return None, f"[{label}]\n{logs}"

    return contigs, f"[{label}]\n{logs}"


def run_assembly_pipeline(
    sample_id,
    r1,
    r2=None,
    outdir=None,
    pre_trimmed: bool = False,
    keep_spades_output: bool = False,
):
    """
    De novo assemble reads with SPAdes. Assembly-only mode with
    isolate optimisation and container-aware resource allocation.
    Quality statistics are computed afterwards.

    Args:
        sample_id: identifier used to name the output directory.
        r1, r2: trimmed read paths (r2 None for single-end).
        outdir: output directory; defaults to RESULTS_DIR/assembly/sample_id.
        pre_trimmed: retained for call-site compatibility.
        keep_spades_output: if True, keeps spades_output/ and its files.

    Returns:
        AssemblyResult(contigs, logs).
        contigs is None if SPAdes failed or was OOM-killed.
    """
    if outdir is None:
        outdir = RESULTS_DIR / "assembly" / sample_id
    else:
        outdir = Path(outdir)

    logger.info(
        "Assembly | %s | start pe=%s",
        sample_id,
        r2 is not None,
    )

    contigs, logs = _run_spades(r1, r2, outdir)

    if contigs and contigs.exists():
        canonical = outdir / f"{sample_id}.fasta"

        if contigs != canonical:
            shutil.copy2(contigs, canonical)
            tool_dir = contigs.parent
            contigs = canonical

            if not keep_spades_output:
                shutil.rmtree(tool_dir, ignore_errors=True)

    if contigs and contigs.exists():
        logger.info(
            "Assembly | %s | done contigs=%s",
            sample_id,
            contigs.name,
        )
    else:
        logger.error(
            "Assembly | %s | failed — no contigs produced",
            sample_id,
        )

    return AssemblyResult(
        contigs=contigs,
        logs=logs,
    )


ASSEMBLY_QC_THRESHOLDS = {
    "min_genome_fraction": 90.0,
    "min_gc": 50.0,
    "max_gc": 56.0,
    "pass": {
        "max_contigs": 150,
        "min_n50": 30_000,
        "min_total_len": 2_000_000,
        "max_total_len": 2_200_000,
    },
    "caution": {
        "max_contigs": 180,
        "min_n50": 20_000,
        "min_total_len": 1_800_000,
        "max_total_len": 2_200_000,
    },
}


def _cdc_tier_hits(stats: dict, tier: dict) -> tuple[int, int, list[str]]:
    hits, total, flags = 0, 0, []

    nc = stats.get("n_contigs")
    if nc is not None:
        total += 1
        ok = nc <= tier["max_contigs"]
        hits += ok
        flags.append(
            f"Contigs {nc} ({'≤' if ok else '>'} {tier['max_contigs']})"
        )

    n50 = stats.get("n50")
    if n50 is not None:
        total += 1
        ok = n50 > tier["min_n50"]
        hits += ok
        flags.append(
            f"N50 {n50/1e3:.1f} kb "
            f"({'>' if ok else '≤'} {tier['min_n50']/1e3:.0f} kb)"
        )

    tl = stats.get("total_len")
    if tl is not None:
        total += 1
        ok = tier["min_total_len"] <= tl <= tier["max_total_len"]
        hits += ok
        flags.append(
            f"Length {tl/1e6:.2f} Mb "
            f"({'within' if ok else 'outside'} "
            f"{tier['min_total_len']/1e6:.1f}–"
            f"{tier['max_total_len']/1e6:.1f} Mb)"
        )

    return hits, total, flags


def evaluate_assembly_qc(stats: dict) -> dict:
    t = ASSEMBLY_QC_THRESHOLDS
    hard_fail_flags: list[str] = []
    hard_fail = False

    gf = stats.get("completeness")
    if gf is not None and gf < t["min_genome_fraction"]:
        hard_fail_flags.append(
            f"Core genes {gf:.1f}% < {t['min_genome_fraction']}%"
        )
        hard_fail = True

    gc = stats.get("gc_pct")
    if gc is not None and not (t["min_gc"] <= gc <= t["max_gc"]):
        hard_fail_flags.append(
            f"GC% {gc:.1f}% outside expected range "
            f"{t['min_gc']}–{t['max_gc']}%"
        )
        hard_fail = True

    if hard_fail:
        return {
            "status": "fail",
            "flags": hard_fail_flags,
        }

    pass_hits, pass_total, pass_flags = _cdc_tier_hits(
        stats,
        t["pass"],
    )

    if pass_total and pass_hits >= min(2, pass_total):
        return {
            "status": "pass",
            "flags": pass_flags,
        }

    caution_hits, caution_total, caution_flags = _cdc_tier_hits(
        stats,
        t["caution"],
    )

    if caution_total and caution_hits >= min(2, caution_total):
        return {
            "status": "caution",
            "flags": caution_flags,
        }

    return {
        "status": "fail",
        "flags": caution_flags or pass_flags,
    }


def parse_contigs_stats(path: str | Path) -> AssemblyStats | None:
    path = Path(path) if path else None
    if not path or not path.exists():
        return None

    lengths = []
    gc_count = 0
    n_count = 0
    total_bases = 0

    for record in SeqIO.parse(str(path), "fasta"):
        s = str(record.seq).upper()
        lengths.append(len(s))
        gc_count += s.count("G") + s.count("C")
        n_count += s.count("N")
        total_bases += len(s)

    if not lengths:
        return None

    lengths.sort(reverse=True)

    def _nx_lx(fraction: float) -> tuple[int, int]:
        cum = 0

        for i, length in enumerate(lengths, start=1):
            cum += length

            if cum >= total_bases * fraction:
                return length, i

        return 0, 0

    n50, l50 = _nx_lx(0.5)
    n90, l90 = _nx_lx(0.9)

    auN = (
        sum(length ** 2 for length in lengths) / total_bases
        if total_bases
        else 0.0
    )

    _stats_d: dict = {
        "n_contigs": len(lengths),
        "total_len": total_bases,
        "n50": n50,
        "n90": n90,
        "l50": l50,
        "l90": l90,
        "auN": auN,
        "largest": lengths[0],
        "gc_pct": (
            gc_count / total_bases * 100
            if total_bases
            else 0.0
        ),
        "contigs_500": sum(
            1 for length in lengths
            if length >= 500
        ),
        "completeness": 0.0,
        "avg_contig_len": (
            total_bases / len(lengths)
            if lengths
            else 0.0
        ),
        "n_per_100kbp": (
            n_count / total_bases * 1e5
            if total_bases
            else 0.0
        ),
    }

    try:
        from backend.phylogeny.cgmlst import (
            cached_core_genome_gene_count
        )
        core = cached_core_genome_gene_count(str(path))
    except Exception:
        core = None

    if core:
        _stats_d["completeness"] = core["core_completeness_pct"]
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
        avg_contig_len=_stats_d.get("avg_contig_len"),
        n_per_100kbp=_stats_d.get("n_per_100kbp"),
        n90=_stats_d.get("n90"),
        l50=_stats_d.get("l50"),
        l90=_stats_d.get("l90"),
        auN=_stats_d.get("auN"),
    )