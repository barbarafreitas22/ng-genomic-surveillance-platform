from pathlib import Path

from .paths import PROJECT_ROOT, NGSTAR_PROFILES_CACHE
from .pubmlst_typing import query_locus_local, load_profiles, lookup_st

ALLELE_DIR = PROJECT_ROOT / "data" / "genes" / "ngstar_alleles"

SCHEME_ID = 67

NGSTAR_GENES = ["penA", "mtrR", "porB", "ponA", "gyrA", "parC", "23SrRNA"]
_PUBMLST_COLUMNS = ["NG_penA", "'mtrR", "NG_porB", "NG_ponA", "NG_gyrA", "NG_parC", "NG_23S"]

_PROFILES_CACHE = NGSTAR_PROFILES_CACHE


def run_ngstar(contigs: Path, minimap2: str = "minimap2") -> dict:
    """
    Run NG-STAR typing on an assembly.
    BLASTs the whole assembly directly against local allele FASTAs.
    """
    alleles = {gene: query_locus_local(gene, contigs, ALLELE_DIR) for gene in NGSTAR_GENES}

    def _valid_allele(v: str) -> bool:
        try:
            float(v)
            return True
        except ValueError:
            return False

    all_valid = all(_valid_allele(v) for v in alleles.values())
    if all_valid:
        profiles = load_profiles(_PROFILES_CACHE, SCHEME_ID, _PUBMLST_COLUMNS)
        st = lookup_st(alleles, NGSTAR_GENES, profiles)
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
