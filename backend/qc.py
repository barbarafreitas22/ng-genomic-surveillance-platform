import json
import logging
import os
import subprocess
from pathlib import Path
import configparser
import zipfile

logger = logging.getLogger(__name__)

_CPU_THREADS = min(os.cpu_count() or 4, 8)

from .paths import PROJECT_ROOT, CONFIG_PATH, RESULTS_DIR, WHOF_REF

config = configparser.ConfigParser()
config.read(CONFIG_PATH)
paths = config["PATHS"]

FASTQC_PATH = paths.get("FASTQC")
FASTP_BIN = paths.get("FASTP", "fastp")
MULTIQC_PATH = paths.get("MULTIQC")
MINIMAP2_BIN = paths.get("MINIMAP2", "minimap2")
KRAKEN2_BIN = paths.get("KRAKEN2", "kraken2")
_kraken2_db_raw = paths.get("KRAKEN2_DB", "")
KRAKEN2_DB = (
    PROJECT_ROOT / _kraken2_db_raw
    if _kraken2_db_raw and not Path(_kraken2_db_raw).is_absolute()
    else Path(_kraken2_db_raw)
)

UPLOAD_DIR = Path(paths.get("UPLOAD_DIR"))


# Neisseria gonorrhoeae NCBI taxonomy ID
_NG_TAXID = "485"
_NG_GENUS_TAXID = "482"


def run_cmd(cmd, step):
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"{step} failed:\n{proc.stderr}")
    return proc.stdout



def run_fastqc(files, outdir, fastqc_path):
    outdir.mkdir(parents=True, exist_ok=True)
    n = len(files)
    run_cmd(
        [fastqc_path, *[str(f) for f in files], "-t", str(n), "-o", str(outdir)],
        "FastQC",
    )


def parse_fastqc(zip_file):
    metrics = {}

    with zipfile.ZipFile(zip_file) as z:
        fname = [f for f in z.namelist() if "fastqc_data.txt" in f][0]
        with z.open(fname) as f:
            for line in f:
                line = line.decode()

                if "Total Sequences" in line:
                    metrics["total_reads"] = int(line.split("\t")[1])

                elif "%GC" in line:
                    metrics["gc"] = int(line.split("\t")[1])

    return metrics


def run_fastp(r1, r2, outdir, fastp_bin="fastp"):
    outdir.mkdir(parents=True, exist_ok=True)
    out_r1 = outdir / "R1.trimmed.fastq.gz"
    out_r2 = outdir / "R2.trimmed.fastq.gz" if r2 else None

    cmd = [
        fastp_bin,
        "--in1", str(r1), "--out1", str(out_r1),
        "--thread", str(min(2, _CPU_THREADS)),
        "--cut_right", "--cut_window_size", "4", "--cut_mean_quality", "20",
        "--length_required", "50",
        "--json", str(outdir / "fastp.json"),
        "--html", str(outdir / "fastp.html"),
    ]

    if r2:
        cmd += ["--in2", str(r2), "--out2", str(out_r2), "--detect_adapter_for_pe"]

    run_cmd(cmd, "fastp")
    return out_r1, out_r2


# MultiQC

def run_multiqc(input_dir, outdir, multiqc_path):
    outdir.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        [multiqc_path, str(input_dir), "-o", str(outdir), "--force"],
        capture_output=True, text=True
    )
    report = outdir / "multiqc_report.html"
    if not report.exists():
        raise RuntimeError(f"MultiQC failed:\n{proc.stderr}")
    return report

def run_kraken2(r1: Path, r2, outdir: Path) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    report = outdir / "kraken2_report.txt"
    output = outdir / "kraken2_output.txt"

    cmd = [
        KRAKEN2_BIN,
        "--db", str(KRAKEN2_DB),
        "--report", str(report),
        "--output", str(output),
        "--threads", str(_CPU_THREADS),
    ]

    if str(r1).endswith(".gz"):
        cmd.append("--gzip-compressed")

    if r2 is not None:
        cmd += ["--paired", str(r1), str(r2)]
    else:
        cmd.append(str(r1))

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Kraken2 failed:\n{proc.stderr}")

    output.unlink(missing_ok=True)
    return report


def parse_kraken2_report(report_path: Path) -> dict:
    """
    Parse Kraken2 report and return species confirmation metrics.
    Columns: pct  reads_clade  reads_direct  rank  taxid  name
    """
    ng_pct        = 0.0
    neisseria_pct = 0.0
    unclassified  = 0.0

    with open(report_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            try:
                pct   = float(parts[0])
                rank  = parts[3].strip()
                taxid = parts[4].strip()
                name  = parts[5].strip()
            except ValueError:
                continue

            if rank == "U":
                unclassified = pct
            if taxid == _NG_TAXID:
                ng_pct = pct
            if taxid == _NG_GENUS_TAXID:
                neisseria_pct = pct

    contamination = ng_pct < 80.0
    if ng_pct >= 95.0:
        species_status = "confirmed"
    elif ng_pct >= 80.0:
        species_status = "likely"
    else:
        species_status = "contaminated"

    return {
        "ng_pct":          round(ng_pct, 2),
        "neisseria_pct":   round(neisseria_pct, 2),
        "unclassified_pct": round(unclassified, 2),
        "species_status":  species_status,
        "contamination":   contamination,
    }


def run_coverage_analysis(r1: Path, r2, ref_genome: Path, outdir: Path) -> dict:
    """Map trimmed reads to reference genome; return mean depth and % genome ≥10×."""
    outdir.mkdir(parents=True, exist_ok=True)
    sorted_bam = outdir / "coverage.sorted.bam"

    mm2_cmd = [MINIMAP2_BIN, "-ax", "sr", "--secondary=no", "-t", str(_CPU_THREADS), str(ref_genome)]
    if r2 is not None and Path(str(r2)).exists():
        mm2_cmd += [str(r1), str(r2)]
    else:
        mm2_cmd.append(str(r1))

    mm2 = subprocess.run(mm2_cmd, capture_output=True)
    if mm2.returncode != 0:
        raise RuntimeError(f"minimap2 coverage failed:\n{mm2.stderr.decode()[:500]}")

    view = subprocess.run(
        ["samtools", "view", "-bS", "-F", "4"],
        input=mm2.stdout, capture_output=True,
    )
    sort = subprocess.run(
        ["samtools", "sort", "-o", str(sorted_bam)],
        input=view.stdout, capture_output=True,
    )
    if sort.returncode != 0:
        raise RuntimeError("samtools sort failed")

    depth_proc = subprocess.run(
        ["samtools", "depth", "-a", str(sorted_bam)],
        capture_output=True, text=True,
    )

    depths = []
    for line in depth_proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3:
            try:
                depths.append(int(parts[2]))
            except ValueError:
                pass

    if not depths:
        return {"mean_coverage": 0.0, "pct_breadth_10x": 0.0}

    mean_cov = sum(depths) / len(depths)
    pct_10x  = sum(1 for d in depths if d >= 10) / len(depths) * 100

    return {
        "mean_coverage":   round(mean_cov, 1),
        "pct_breadth_10x": round(pct_10x, 1),
    }


def _run_kraken2_safe(r1, r2, outdir) -> tuple:
    if not KRAKEN2_DB.exists():
        return {"species_status": "db_not_found"}, None
    try:
        report_path = run_kraken2(r1, r2, outdir)
        return parse_kraken2_report(report_path), report_path
    except Exception as e:
        return {"error": str(e), "species_status": "error"}, None


def run_qc_pipeline_fasta(sample_id, fasta_path, output_dir=None):
    """Kraken2 species confirmation for an already-assembled genome. No trimming."""
    fasta_path = Path(fasta_path)
    base = Path(output_dir) if output_dir else RESULTS_DIR / sample_id
    base.mkdir(parents=True, exist_ok=True)

    logger.info("QC | %s | start (FASTA)", sample_id)

    kraken2_metrics, kraken2_report_path = _run_kraken2_safe(fasta_path, None, base / "kraken2")
    logger.info("QC | %s | kraken2 species_status=%s ng_pct=%.1f%%",
                sample_id,
                kraken2_metrics.get("species_status", "?"),
                kraken2_metrics.get("ng_pct", 0.0))

    return {
        "metrics":        {},
        "kraken2":        kraken2_metrics,
        "kraken2_report": str(kraken2_report_path) if kraken2_report_path else None,
        "trimmed_r1":     None,
        "trimmed_r2":     None,
        "multiqc":        None,
        "job_dir":        str(base),
        "input_type":     "fasta",
    }


def run_qc_pipeline(sample_id, file_list, output_dir=None):
    """
    QC pipeline: FastQC → fastp → Kraken2 → MultiQC
    file_list: list of Path objects (1 file = SE, 2 files = PE)
    """
    file_list = [Path(f) for f in (file_list if isinstance(file_list, list) else [file_list])]
    r1 = file_list[0]
    r2 = file_list[1] if len(file_list) >= 2 else None

    base = Path(output_dir) if output_dir else RESULTS_DIR / sample_id
    job_dir = base / Path(r1).stem
    job_dir.mkdir(parents=True, exist_ok=True)

    run_fastqc(
        [r1] if r2 is None else [r1, r2],
        job_dir / "fastqc",
        FASTQC_PATH,
    )

    fastqc_zips = list((job_dir / "fastqc").glob("*.zip"))
    fastqc_metrics = parse_fastqc(fastqc_zips[0]) if fastqc_zips else {}

    r1_trim, r2_trim = run_fastp(r1, r2, job_dir / "trim", FASTP_BIN)

    kraken2_metrics, kraken2_report_path = _run_kraken2_safe(r1_trim, r2_trim, job_dir / "kraken2")

    coverage_metrics: dict = {}
    if WHOF_REF.exists():
        try:
            coverage_metrics = run_coverage_analysis(r1_trim, r2_trim, WHOF_REF, job_dir / "coverage")
        except Exception as e:
            logger.warning("QC | %s | coverage failed: %s", sample_id, e)

    multiqc_report = run_multiqc(
        job_dir,
        job_dir / "multiqc",
        MULTIQC_PATH,
    )

    return {
        "metrics":    {**fastqc_metrics, **kraken2_metrics, **coverage_metrics},
        "kraken2":    kraken2_metrics,
        "kraken2_report": str(kraken2_report_path) if kraken2_report_path else None,
        "trimmed_r1": str(r1_trim),
        "trimmed_r2": str(r2_trim) if r2_trim else None,
        "multiqc":    str(multiqc_report),
        "job_dir":    str(job_dir),
    }
