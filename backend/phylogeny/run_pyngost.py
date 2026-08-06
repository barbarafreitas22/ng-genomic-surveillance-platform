import csv
import glob
import os
import subprocess
import sys
from pathlib import Path

from ..paths import PYNGOST_DB_DIR as _PYNGOST_DB_DIR
PYNGOST_DB_DIR = str(_PYNGOST_DB_DIR)


def _script_path() -> str:
    for p in sys.path:
        candidate = os.path.join(p, "pyngoST", "pyngoST.py")
        if os.path.exists(candidate):
            return candidate
    raise RuntimeError("pyngoST not found — install with: pip install pyngoST")


def is_pyngost_ready() -> bool:
    try:
        _script_path()
    except RuntimeError:
        return False
    return os.path.isdir(PYNGOST_DB_DIR) and bool(
        glob.glob(os.path.join(PYNGOST_DB_DIR, "*.fas"))
    )


def download_pyngost_db() -> str:
    """
    Download MLST / NG-STAR / NG-MAST allele databases from PubMLST and ngstar.canada.ca.
    Only needs to run once; results are cached in PYNGOST_DB_DIR.
    """
    os.makedirs(PYNGOST_DB_DIR, exist_ok=True)
    script = _script_path()
    proc = subprocess.run(
        [sys.executable, script, "-d", "-n", PYNGOST_DB_DIR],
        cwd=os.path.dirname(script),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if not glob.glob(os.path.join(PYNGOST_DB_DIR, "*.fas")):
        raise RuntimeError(f"Database download failed:\n{proc.stdout[-2000:]}")
    return PYNGOST_DB_DIR


def run_typing(genome_paths: list[str], run_dir: str) -> dict[str, dict]:
    """
    Run MLST + NG-STAR + NG-MAST multi-scheme typing on assembled FASTA files.

    Returns:
        {sample_stem: {mlst_st, ngstar_st, ngstar_cc, mosaic_pena, ngmast_st, ngmast_genogroup}}

    Keys match the FASTA file stems passed in genome_paths.
    """
    if not is_pyngost_ready():
        raise RuntimeError(
            f"pyngoST allele database not found at {PYNGOST_DB_DIR}. "
            "Download it first from the Phylogenetics page."
        )

    out_dir = os.path.join(run_dir, "pyngost_out")
    os.makedirs(out_dir, exist_ok=True)
    out_file = "typing_results.tsv"
    script = _script_path()

    cmd = [
        sys.executable, script,
        "-i", *genome_paths,
        "-s", "NG-STAR,MLST,NG-MAST",
        "-g",       # NG-MAST genogroups
        "-c",       # NG-STAR clonal complexes
        "-m",       # mosaic/semimosaic penA detection
        "-p", PYNGOST_DB_DIR,
        "-q", out_dir,
        "-o", out_file,
        "-t", "4",
    ]

    proc = subprocess.run(
        cmd,
        cwd=os.path.dirname(script),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    result_path = os.path.join(out_dir, out_file)
    if not os.path.exists(result_path):
        raise RuntimeError(f"pyngoST produced no output:\n{proc.stdout[-2000:]}")

    results: dict[str, dict] = {}
    with open(result_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            strain = row.get("strain", "")
            stem = Path(strain).stem or strain
            results[stem] = {
                "mlst_st":          row.get("MLST", "-"),
                "ngstar_st":        row.get("NG-STAR", "-"),
                "ngstar_cc":        row.get("NG-STAR_CC", "-"),
                "mosaic_pena":      row.get("penA_mosaic_type", "-"),
                "ngmast_st":        row.get("NG-MAST", "-"),
                "ngmast_genogroup": row.get("Genogroup", "-"),
            }
    return results
