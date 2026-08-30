import json
import logging
import os
import subprocess
import tempfile
from pathlib import Path
import configparser

logger = logging.getLogger(__name__)

_CPU_THREADS = min(os.cpu_count() or 4, 8)
_FASTP_THREADS = min(4, _CPU_THREADS)
_SAMTOOLS_THREADS = min(4, _CPU_THREADS)

from .paths import PROJECT_ROOT, CONFIG_PATH, RESULTS_DIR, WHOF_REF


config = configparser.ConfigParser()
config.read(CONFIG_PATH)
paths = config["PATHS"]

FASTP_BIN = paths.get("FASTP", "fastp")
MINIMAP2_BIN = paths.get("MINIMAP2", "minimap2")
KRAKEN2_BIN = paths.get("KRAKEN2", "kraken2")

_kraken2_db_raw = paths.get("KRAKEN2_DB", "")
KRAKEN2_DB = (
    PROJECT_ROOT / _kraken2_db_raw
    if _kraken2_db_raw and not Path(_kraken2_db_raw).is_absolute()
    else Path(_kraken2_db_raw)
)

UPLOAD_DIR = Path(paths.get("UPLOAD_DIR"))


_NG_TAXID = "485"
_NG_GENUS_TAXID = "482"


def run_cmd(cmd, step):
    proc = subprocess.run(cmd, capture_output=True, text=True)

    if proc.returncode != 0:
        if proc.returncode in (-9, 137):
            raise RuntimeError(
                f"{step} killed by OOM (rc={proc.returncode}) — "
                f"cgroup memory limit reached.\n{proc.stderr}"
            )

        raise RuntimeError(
            f"{step} failed (rc={proc.returncode}):\n{proc.stderr}"
        )

    return proc.stdout


def parse_fastp_json(json_path):
    """
    Extract pre-trimming read count and GC content from fastp's
    before_filtering summary.
    """
    metrics = {}

    with open(json_path) as f:
        data = json.load(f)

    before = data.get("summary", {}).get("before_filtering", {})

    if "total_reads" in before:
        metrics["total_reads"] = before["total_reads"]

    if "gc_content" in before:
        metrics["gc"] = round(before["gc_content"] * 100)

    return metrics


def run_fastp(r1, r2, outdir, fastp_bin="fastp"):
    outdir.mkdir(parents=True, exist_ok=True)

    out_r1 = outdir / "R1.trimmed.fastq.gz"
    out_r2 = outdir / "R2.trimmed.fastq.gz" if r2 else None

    cmd = [
        fastp_bin,
        "--in1", str(r1),
        "--out1", str(out_r1),

        "--thread", str(_FASTP_THREADS),

        "--cut_right",
        "--cut_window_size", "4",
        "--cut_mean_quality", "20",
        "--length_required", "50",

        "--json", str(outdir / "fastp.json"),
        "--html", str(outdir / "fastp.html"),
    ]

    if r2:
        cmd += [
            "--in2", str(r2),
            "--out2", str(out_r2),
            "--detect_adapter_for_pe",
        ]

    run_cmd(cmd, "fastp")

    return out_r1, out_r2



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

    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
    )

    if proc.returncode != 0:
        raise RuntimeError(f"Kraken2 failed:\n{proc.stderr}")

    output.unlink(missing_ok=True)

    return report


def parse_kraken2_report(report_path: Path) -> dict:
    """
    Parse a Kraken2 report and classify the sample by species purity.

    Taxonomy:
        485 = N. gonorrhoeae
        482 = Neisseria genus

    Species status:
        >=95% -> confirmed
        >=80% -> likely
        <80%  -> contaminated
    """
    ng_pct = 0.0
    neisseria_pct = 0.0
    unclassified = 0.0

    with open(report_path) as f:
        for line in f:
            parts = line.strip().split("\t")

            if len(parts) < 6:
                continue

            try:
                pct = float(parts[0])
                rank = parts[3].strip()
                taxid = parts[4].strip()
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
        "ng_pct": round(ng_pct, 2),
        "neisseria_pct": round(neisseria_pct, 2),
        "unclassified_pct": round(unclassified, 2),
        "species_status": species_status,
        "contamination": contamination,
    }



def run_coverage_analysis(
    r1: Path,
    r2,
    ref_genome: Path,
    outdir: Path,
) -> dict:
    """
    Map trimmed reads to a reference genome and calculate coverage metrics.

    Pipeline:
        minimap2 -> samtools view -> samtools sort -> samtools depth

    The alignment steps are connected through pipes so that the full SAM/BAM
    stream is not loaded into Python memory.

    Returns:
        mean_coverage:
            Mean sequencing depth across the complete reference.

        pct_breadth_10x:
            Percentage of reference positions covered at >=10x.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    sorted_bam = outdir / "coverage.sorted.bam"

    mm2_cmd = [
        MINIMAP2_BIN,
        "-ax", "sr",
        "--secondary=no",
        "-t", str(_CPU_THREADS),
        str(ref_genome),
    ]

    if r2 is not None and Path(str(r2)).exists():
        mm2_cmd += [str(r1), str(r2)]
    else:
        mm2_cmd.append(str(r1))

    view_cmd = [
        "samtools",
        "view",
        "-@",
        str(_SAMTOOLS_THREADS),
        "-b",
        "-F",
        "4",
    ]

    sort_cmd = [
        "samtools",
        "sort",
        "-@",
        str(_SAMTOOLS_THREADS),
        "-o",
        str(sorted_bam),
    ]

    with tempfile.TemporaryFile(mode="w+b") as log_file:

        mm2 = subprocess.Popen(
            mm2_cmd,
            stdout=subprocess.PIPE,
            stderr=log_file,
        )

        view = subprocess.Popen(
            view_cmd,
            stdin=mm2.stdout,
            stdout=subprocess.PIPE,
            stderr=log_file,
        )

        if mm2.stdout is not None:
            mm2.stdout.close()

        sort = subprocess.Popen(
            sort_cmd,
            stdin=view.stdout,
            stdout=subprocess.DEVNULL,
            stderr=log_file,
        )

        if view.stdout is not None:
            view.stdout.close()

        sort_rc = sort.wait()
        view_rc = view.wait()
        mm2_rc = mm2.wait()

        log_file.seek(0)
        stderr_text = log_file.read().decode(
            errors="replace"
        ).strip()

    if mm2_rc != 0:
        if mm2_rc in (-9, 137):
            raise RuntimeError(
                f"minimap2 coverage killed by OOM "
                f"(rc={mm2_rc}).\n{stderr_text}"
            )
        raise RuntimeError(
            f"minimap2 coverage failed (rc={mm2_rc}):\n{stderr_text}"
        )

    if view_rc != 0:
        if view_rc in (-9, 137):
            raise RuntimeError(
                f"samtools view killed by OOM "
                f"(rc={view_rc}).\n{stderr_text}"
            )
        raise RuntimeError(
            f"samtools view failed (rc={view_rc}):\n{stderr_text}"
        )

    if sort_rc != 0:
        if sort_rc in (-9, 137):
            raise RuntimeError(
                f"samtools sort killed by OOM "
                f"(rc={sort_rc}).\n{stderr_text}"
            )
        raise RuntimeError(
            f"samtools sort failed (rc={sort_rc}):\n{stderr_text}"
        )

    depth_cmd = [
        "samtools",
        "depth",
        "-a",
        str(sorted_bam),
    ]

    total_positions = 0
    covered_10x = 0
    depth_sum = 0

    with subprocess.Popen(
        depth_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    ) as depth_proc:

        if depth_proc.stdout is not None:
            for line in depth_proc.stdout:
                parts = line.split("\t")

                if len(parts) < 3:
                    continue

                try:
                    depth = int(parts[2])
                except ValueError:
                    continue

                total_positions += 1
                depth_sum += depth

                if depth >= 10:
                    covered_10x += 1

        depth_stderr = (
            depth_proc.stderr.read()
            if depth_proc.stderr is not None
            else ""
        )

        depth_rc = depth_proc.wait()

    if depth_rc != 0:
        if depth_rc in (-9, 137):
            raise RuntimeError(
                f"samtools depth killed by OOM "
                f"(rc={depth_rc}).\n{depth_stderr}"
            )

        raise RuntimeError(
            f"samtools depth failed (rc={depth_rc}):\n"
            f"{depth_stderr}"
        )

    if total_positions == 0:
        return {
            "mean_coverage": 0.0,
            "pct_breadth_10x": 0.0,
        }

    mean_cov = depth_sum / total_positions
    pct_10x = covered_10x / total_positions * 100

    return {
        "mean_coverage": round(mean_cov, 1),
        "pct_breadth_10x": round(pct_10x, 1),
    }


def _run_kraken2_safe(r1, r2, outdir) -> tuple:
    if not KRAKEN2_DB.exists():
        return {"species_status": "db_not_found"}, None

    try:
        report_path = run_kraken2(r1, r2, outdir)
        return parse_kraken2_report(report_path), report_path

    except Exception as e:
        return {
            "error": str(e),
            "species_status": "error",
        }, None


def run_qc_pipeline_fasta(
    sample_id,
    fasta_path,
    output_dir=None,
):
    """
    Species-confirmation QC for an already-assembled genome.

    No read trimming or coverage analysis is performed because no raw reads
    are available.
    """
    fasta_path = Path(fasta_path)

    base = (
        Path(output_dir)
        if output_dir
        else RESULTS_DIR / sample_id
    )

    base.mkdir(parents=True, exist_ok=True)

    logger.info(
        "QC | %s | start (FASTA)",
        sample_id,
    )

    kraken2_metrics, kraken2_report_path = _run_kraken2_safe(
        fasta_path,
        None,
        base / "kraken2",
    )

    logger.info(
        "QC | %s | kraken2 species_status=%s ng_pct=%.1f%%",
        sample_id,
        kraken2_metrics.get("species_status", "?"),
        kraken2_metrics.get("ng_pct", 0.0),
    )

    return {
        "metrics": {},
        "kraken2": kraken2_metrics,
        "kraken2_report": (
            str(kraken2_report_path)
            if kraken2_report_path
            else None
        ),
        "trimmed_r1": None,
        "trimmed_r2": None,
        "job_dir": str(base),
        "input_type": "fasta",
    }


def run_qc_pipeline(
    sample_id,
    file_list,
    output_dir=None,
):
    """
    Raw-read QC pipeline:

        fastp
          -> Kraken2
          -> minimap2/samtools coverage

    Entry point for Module 1 when raw FASTQ reads are uploaded.

    Pre-trimming read count and GC content are obtained from fastp's own
    report. A separate FastQC pass is intentionally omitted because it
    would duplicate these basic read-level metrics while adding another
    complete scan of the FASTQ files.

    Args:
        sample_id:
            Sample identifier.

        file_list:
            One Path for single-end or two Paths [R1, R2] for paired-end.

        output_dir:
            Base output directory. Defaults to RESULTS_DIR/sample_id.

    Returns:
        Dictionary containing:
            - merged QC metrics
            - trimmed FASTQ paths
            - Kraken2 report path
            - job directory
    """
    file_list = [
        Path(f)
        for f in (
            file_list
            if isinstance(file_list, list)
            else [file_list]
        )
    ]

    r1 = file_list[0]
    r2 = file_list[1] if len(file_list) >= 2 else None

    base = (
        Path(output_dir)
        if output_dir
        else RESULTS_DIR / sample_id
    )

    job_dir = base / Path(r1).stem
    job_dir.mkdir(parents=True, exist_ok=True)


    trim_dir = job_dir / "trim"

    r1_trim, r2_trim = run_fastp(
        r1,
        r2,
        trim_dir,
        FASTP_BIN,
    )

    fastp_json = trim_dir / "fastp.json"

    pretrim_metrics = (
        parse_fastp_json(fastp_json)
        if fastp_json.exists()
        else {}
    )


    kraken2_metrics, kraken2_report_path = _run_kraken2_safe(
        r1_trim,
        r2_trim,
        job_dir / "kraken2",
    )


    coverage_metrics: dict = {}

    if WHOF_REF.exists():
        try:
            coverage_metrics = run_coverage_analysis(
                r1_trim,
                r2_trim,
                WHOF_REF,
                job_dir / "coverage",
            )

        except Exception as e:
            logger.warning(
                "QC | %s | coverage failed: %s",
                sample_id,
                e,
            )



    return {
        "metrics": {
            **pretrim_metrics,
            **kraken2_metrics,
            **coverage_metrics,
        },
        "kraken2": kraken2_metrics,
        "kraken2_report": (
            str(kraken2_report_path)
            if kraken2_report_path
            else None
        ),
        "trimmed_r1": str(r1_trim),
        "trimmed_r2": (
            str(r2_trim)
            if r2_trim
            else None
        ),
        "job_dir": str(job_dir),
    }
