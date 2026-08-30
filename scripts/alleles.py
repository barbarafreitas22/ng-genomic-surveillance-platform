"""
Download all allele FASTA sequences from PubMLST for NG-STAR and MLST.
"""
import logging
import sys
import time
from pathlib import Path

import requests

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s  %(message)s",
    stream=sys.stdout,
)
log = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
API_BASE = "https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef"


NGSTAR_GENES = {
    "mtrR":    "'mtrR",
    "porB":    "NG_porB",
    "ponA":    "NG_ponA",
    "gyrA":    "NG_gyrA",
    "parC":    "NG_parC",
    "23SrRNA": "NG_23S",
}

MLST_GENES = {
    "abcZ": "abcZ",
    "adk":  "adk",
    "aroE": "aroE",
    "fumC": "fumC",
    "gdh":  "gdh",
    "pdhC": "pdhC",
    "pgm":  "pgm",
}

NGMAST_GENES = {
    "porB": "NG-MAST_porB",
    "tbpB": "NG-MAST_tbpB",
}


def download_alleles(locus: str, out_path: Path, retries: int = 3) -> bool:
    url = f"{API_BASE}/loci/{locus}/alleles_fasta"
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, timeout=120)
            if r.status_code == 200 and r.text.strip().startswith(">"):
                out_path.write_text(r.text)
                n = r.text.count(">")
                log.info("%-20s → %s  (%d alleles)", locus, out_path.name, n)
                return True
            log.warning("%s — HTTP %s (attempt %d)", locus, r.status_code, attempt)
        except requests.RequestException as e:
            log.warning("%s — %s (attempt %d)", locus, e, attempt)
        if attempt < retries:
            time.sleep(2)
    return False


def main():
    ngstar_dir = PROJECT_ROOT / "data" / "genes" / "ngstar_alleles"
    mlst_dir   = PROJECT_ROOT / "data" / "genes" / "mlst_alleles"
    ngmast_dir = PROJECT_ROOT / "data" / "genes" / "ngmast_alleles"
    ngstar_dir.mkdir(parents=True, exist_ok=True)
    mlst_dir.mkdir(parents=True, exist_ok=True)
    ngmast_dir.mkdir(parents=True, exist_ok=True)

    log.info("=== NG-STAR alleles ===")
    ok = True
    for gene, locus in NGSTAR_GENES.items():
        ok &= download_alleles(locus, ngstar_dir / f"{gene}.fasta")
        time.sleep(0.5)

    log.info("=== MLST alleles ===")
    for gene, locus in MLST_GENES.items():
        ok &= download_alleles(locus, mlst_dir / f"{gene}.fasta")
        time.sleep(0.5)

    log.info("=== NG-MAST alleles ===")
    for gene, locus in NGMAST_GENES.items():
        ok &= download_alleles(locus, ngmast_dir / f"{gene}.fasta")
        time.sleep(0.5)

    if ok:
        log.info("All downloads completed successfully.")
    else:
        log.error("Downloads failed.")
        sys.exit(1)


if __name__ == "__main__":
    main()
