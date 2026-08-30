from pathlib import Path

from .paths import PROJECT_ROOT, NGMAST_PROFILES_CACHE
from .pubmlst_typing import query_locus_local, load_profiles, lookup_st

ALLELE_DIR = PROJECT_ROOT / "data" / "genes" / "ngmast_alleles"

SCHEME_ID = 71

NGMAST_GENES = ["porB", "tbpB"]
_PUBMLST_COLUMNS = ["NG-MAST_porB", "NG-MAST_tbpB"]

_PROFILES_CACHE = NGMAST_PROFILES_CACHE


def run_ngmast(contigs: Path) -> dict:
    """
    NG-MAST typing against the local PubMLST allele set.
    Args:
        contigs: path to the assembled genome FASTA.

    Returns:
        dict with ST (NG-MAST sequence type), alleles, and an
        incomplete/novel-allele flag.
    """
    alleles = {gene: query_locus_local(gene, contigs, ALLELE_DIR) for gene in NGMAST_GENES}

    def _valid_allele(v: str) -> bool:
        try:
            float(v)
            return True
        except ValueError:
            return False

    all_valid = all(_valid_allele(v) for v in alleles.values())
    if all_valid:
        profiles = load_profiles(_PROFILES_CACHE, SCHEME_ID, _PUBMLST_COLUMNS)
        st = lookup_st(alleles, NGMAST_GENES, profiles)
    else:
        st = "?"

    has_novel = any(v == "new" for v in alleles.values())
    has_error = any(v in ("error", "not_found", "no_ref", "no_db", "?") for v in alleles.values())

    return {
        "ST":         st,
        "alleles":    alleles,
        "novel":      has_novel,
        "incomplete": has_error,
    }
