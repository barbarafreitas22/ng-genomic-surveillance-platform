import subprocess
from pathlib import Path

from .paths import PROJECT_ROOT, MLST_PROFILES_CACHE
from .pubmlst_typing import query_locus_local, load_profiles, lookup_st

ALLELE_DIR = PROJECT_ROOT / "data" / "genes" / "mlst_alleles"

SCHEME_ID  = 1
MLST_GENES = ["abcZ", "adk", "aroE", "fumC", "gdh", "pdhC", "pgm"]

_PROFILES_CACHE = MLST_PROFILES_CACHE


def run_mlst(contigs: Path) -> dict:
    """
    7-gene MLST typing (abcZ, adk, aroE, fumC, gdh, pdhC, pgm) against
    the local PubMLST allele set (pubmlst_typing.query_locus_local),
    falling back to the mlst CLI tool if no local allele DB is present.

    Args:
        contigs: path to the assembled genome FASTA.

    Returns:
        dict with st (sequence type, or "Novel"/"?"), alleles (per-gene
        allele numbers), method ("local" or "cli"), and incomplete/novel
        flags.
    """
    if not ALLELE_DIR.exists() or not any(ALLELE_DIR.glob("*.fasta")):
        try:
            proc = subprocess.run(
                ["mlst", "--scheme", "neisseria", "--quiet", str(contigs)],
                capture_output=True, text=True, timeout=180,
            )
            if proc.returncode == 0:
                lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
                if lines:
                    parts = lines[-1].split("\t")
                    if len(parts) >= 3:
                        st = parts[2].strip()
                        novel = st in ("-", "~") or "~" in st or "?" in st
                        alleles: dict = {}
                        for i, gene in enumerate(MLST_GENES):
                            col = i + 3
                            if col < len(parts):
                                token = parts[col].strip()
                                num = token.split("(")[-1].rstrip(")") if "(" in token else token
                                alleles[gene] = num
                        return {
                            "st":       st if not novel else f"Novel ({st})",
                            "novel":    novel,
                            "alleles":  alleles,
                            "scheme":   "neisseria",
                            "method":   "cli",
                            "incomplete": False,
                        }
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass
        return {"error": "No local MLST alleles — run scripts/download_allele_fastas.py", "st": None, "alleles": {}}

    alleles = {gene: query_locus_local(gene, contigs, ALLELE_DIR) for gene in MLST_GENES}

    all_numeric = all(v.isdigit() for v in alleles.values())
    if all_numeric:
        profiles = load_profiles(_PROFILES_CACHE, SCHEME_ID, MLST_GENES)
        st = lookup_st(alleles, MLST_GENES, profiles)
    else:
        st = "?"
    novel       = any(v == "new" for v in alleles.values())
    incomplete  = any(v in ("error", "not_found", "no_ref", "no_db") for v in alleles.values())

    return {
        "st":        st if all_numeric else ("Novel" if novel else "?"),
        "novel":     novel,
        "alleles":   alleles,
        "scheme":    "neisseria",
        "method":    "local",
        "incomplete": incomplete,
    }
