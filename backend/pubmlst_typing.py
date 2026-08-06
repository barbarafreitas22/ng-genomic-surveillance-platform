import csv
import io
import subprocess
import time
from pathlib import Path
from typing import Optional

import requests

from .paths import PROJECT_ROOT

PUBMLST_DB = "pubmlst_neisseria_seqdef"
API_BASE   = f"https://rest.pubmlst.org/db/{PUBMLST_DB}"
_CACHE_TTL_DAYS = 30


def pubmlst_headers() -> dict:
    key_file = PROJECT_ROOT / "backend" / ".pubmlst_key"
    if key_file.exists():
        key = key_file.read_text().strip()
        if key:
            return {"X-API-Key": key}
    return {}


def query_locus_local(gene: str, assembly: Path, allele_dir: Path) -> str:
    """
    BLAST whole assembly directly against the allele FASTA.
    The coverage is measured against the allele (scov).

    Thresholds:
    >= 99.9 and scov >= 0.95  -> exact allele
    >= 95.0 and scov >= 0.80  -> closest known
      below                  -> new
    """
    alleles_fasta = allele_dir / f"{gene}.fasta"
    if not alleles_fasta.exists():
        return "no_db"

    try:
        proc = subprocess.run(
            [
                "blastn",
                "-query",           str(assembly),
                "-subject",         str(alleles_fasta),
                "-outfmt",          "6 sseqid pident length slen",
                "-perc_identity",   "90",
                "-max_target_seqs", "10",
                "-dust",            "no",
            ],
            capture_output=True, text=True, timeout=120,
        )
    except FileNotFoundError:
        return "error"
    except subprocess.TimeoutExpired:
        return "error"

    if proc.returncode != 0 or not proc.stdout.strip():
        return "new"

    best: Optional[tuple] = None
    best_score = -1.0

    for line in proc.stdout.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        try:
            sseqid = parts[0].split()[0]
            pident = float(parts[1])
            length = int(parts[2])
            slen   = int(parts[3])
        except (ValueError, IndexError):
            continue
        scov  = length / slen if slen > 0 else 0.0
        score = pident * scov
        if score > best_score:
            best_score = score
            best = (sseqid, pident, scov)

    if best is None:
        return "new"

    sseqid, pident, scov = best
    # PubMLST FASTA headers: >abcZ_1, >NG_penA_14, >'mtrR_3, ...
    allele_id = sseqid.rsplit("_", 1)[-1]

    if pident >= 99.9 and scov >= 0.95:
        return allele_id
    if pident >= 95.0 and scov >= 0.80:
        return f"~{allele_id}"
    return "new"


def _parse_profiles_csv(text: str, columns: list[str]) -> dict[tuple, str]:
    profiles: dict[tuple, str] = {}
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    for row in reader:
        st = row.get("ST", "").strip()
        if not st:
            continue
        key = tuple(row.get(col, "?").strip() for col in columns)
        profiles[key] = st
    return profiles


def load_profiles(cache_path: Path, scheme_id: int, columns: list[str]) -> dict[tuple, str]:
    """Load ST profiles for a PubMLST scheme, using a local TSV cache (30-day TTL)."""
    if cache_path.exists():
        age_days = (time.time() - cache_path.stat().st_mtime) / 86400
        if age_days < _CACHE_TTL_DAYS:
            profiles = _parse_profiles_csv(cache_path.read_text(), columns)
            if profiles:
                return profiles

    url = f"{API_BASE}/schemes/{scheme_id}/profiles_csv"
    try:
        resp = requests.get(url, headers=pubmlst_headers(), timeout=60)
        resp.raise_for_status()
        cache_path.write_text(resp.text)
        return _parse_profiles_csv(resp.text, columns)
    except Exception:
        if cache_path.exists():
            return _parse_profiles_csv(cache_path.read_text(), columns)
        return {}


def lookup_st(alleles: dict[str, str], gene_order: list[str], profiles: dict[tuple, str]) -> str:
    if not profiles:
        return "?"
    key = tuple(alleles.get(gene, "?") for gene in gene_order)
    return profiles.get(key, "new")
