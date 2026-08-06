import json
import logging
import os
import re as _re
import subprocess
import threading
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from pathlib import Path
import configparser

from Bio.Seq import Seq

logger = logging.getLogger(__name__)

from .paths import PROJECT_ROOT, CONFIG_PATH, GENE_DB_CROM, GENE_DB_PLASM, GENE_DB_ESSENTIAL

config = configparser.ConfigParser()
config.read(CONFIG_PATH)
paths = config["PATHS"]

MINIMAP2 = paths.get("MINIMAP2", "minimap2")
BLASTN   = paths.get("BLASTN", "blastn")
_FREQ_PANEL_PATH = PROJECT_ROOT / "data" / "synonymous_freq_panel.json"

_FAIDX_LOCK = threading.Lock()

@lru_cache(maxsize=None)
def _samtools_sort_major() -> int:
    r = subprocess.run(["samtools", "--version"], capture_output=True, text=True)
    if r.returncode == 0:
        m = _re.search(r"samtools (\d+)\.", r.stdout)
        if m:
            return int(m.group(1))
    return 0

_RNA_GENES = {"23SrRNA", "16SrRNA"}

_RNA_GENE_POSITIONS: dict[str, frozenset] = {
    "23SrRNA": frozenset({2045, 2597}),
    "16SrRNA": frozenset({1192}),
}

# mtrR_promoter: handled separately via detect_mtrr_promoter_deletion 
_SKIP_CHROM_GENES = {
    "mtrR_promoter",
    "mtrR_promoter_mosaic_1", "mtrR_promoter_mosaic_2", "mtrR_promoter_mosaic_3",
    "mtrD_mosaic_1", "mtrD_mosaic_2", "mtrD_mosaic_3",
}

def _translate_codon(codon: str) -> str:
    try:
        return str(Seq(codon).translate(table=11))
    except Exception:
        return "X"


def _read_fasta_seq(fasta_path: Path) -> str:
    seq = []
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()


def _rc(s: str) -> str:
    return str(Seq(s).reverse_complement())


def _find(probe: str, seq: str) -> bool:
    return probe in seq or _rc(probe) in seq


_AA_POSITION_OFFSET: dict[str, int] = {}


CDC_RULES = {
    "penA_A501V": "ceftriaxone_reduced_susceptibility",
    "penA_A501P": "ceftriaxone_reduced_susceptibility",
    "penA_A501T": "ceftriaxone_reduced_susceptibility",
    "penA_G545S": "ceftriaxone_reduced_susceptibility",
    "penA_N512Y": "ceftriaxone_reduced_susceptibility",
    "penA_P551S": "ceftriaxone_reduced_susceptibility",
    "penA_P551L": "ceftriaxone_reduced_susceptibility",
    "penA_A311V": "reduced_beta_lactam_susceptibility",
    "penA_I312M": "reduced_beta_lactam_susceptibility",
    "penA_V316T": "reduced_beta_lactam_susceptibility",
    "penA_V316P": "reduced_beta_lactam_susceptibility",
    "penA_T483S": "reduced_beta_lactam_susceptibility",
    "penA_F504L": "reduced_beta_lactam_susceptibility",
    "penA_G542S": "reduced_beta_lactam_susceptibility",
    "penA_T534A": "reduced_beta_lactam_susceptibility",
    "penA_A510V": "reduced_beta_lactam_susceptibility",
    "penA_A516G": "reduced_beta_lactam_susceptibility",
    "penA_H541N": "reduced_beta_lactam_susceptibility",
    "penA_P552V": "reduced_beta_lactam_susceptibility",
    "penA_K555Q": "reduced_beta_lactam_susceptibility",
    "penA_I556V": "reduced_beta_lactam_susceptibility",
    "penA_I566V": "reduced_beta_lactam_susceptibility",
    "penA_A574V": "reduced_beta_lactam_susceptibility",
    "penA_ins345": "reduced_beta_lactam_susceptibility",
    "penA_ins346": "reduced_beta_lactam_susceptibility",
    "penA_ins573": "reduced_beta_lactam_susceptibility",

    "23SrRNA_A2045G": "high_level_azithromycin_resistance",
    "23SrRNA_C2597T": "moderate_azithromycin_resistance",

    "16SrRNA_C1192T": "spectinomycin_resistance",

    "parE_G410V": "fluoroquinolone_minor_parE",
    "parE_D420N": "fluoroquinolone_minor_parE",
    "parE_L445H": "fluoroquinolone_minor_parE",
    "parE_T451I": "fluoroquinolone_minor_parE",
    "parE_T451S": "fluoroquinolone_minor_parE",

    "gyrB_G447V": "fluoroquinolone_minor_gyrB",
    "gyrB_D429N": "fluoroquinolone_minor_gyrB",
    "gyrB_A508V": "fluoroquinolone_minor_gyrB",
    "gyrB_K450T": "fluoroquinolone_minor_gyrB",
    "gyrB_E466D": "fluoroquinolone_minor_gyrB",

    "gyrA_A67P": "ciprofloxacin_resistant",
    "gyrA_S91F": "ciprofloxacin_resistant",
    "gyrA_S91Y": "ciprofloxacin_resistant",
    "gyrA_D95G": "ciprofloxacin_resistant",
    "gyrA_D95N": "ciprofloxacin_resistant",
    "gyrA_D95A": "ciprofloxacin_resistant",
    "gyrA_D95E": "ciprofloxacin_resistant",

    "parC_D86N": "ciprofloxacin_resistant",
    "parC_S87R": "ciprofloxacin_resistant",
    "parC_S87W": "ciprofloxacin_resistant",
    "parC_A87G": "ciprofloxacin_resistant",
    "parC_S87N": "ciprofloxacin_intermediate_parC",
    "parC_S87I": "ciprofloxacin_intermediate_parC",
    "parC_S88P": "ciprofloxacin_intermediate_parC",
    "parC_E91K": "ciprofloxacin_intermediate_parC",

    "mtrR_A39T":            "efflux_pump_overexpression",
    "mtrR_G45D":            "efflux_pump_overexpression",
    "mtrR_promoter_del35A": "efflux_pump_overexpression",
    "mtrR_promoter_AtoC":   "efflux_pump_overexpression",
    "mtrR_promoter_mtr120": "efflux_pump_overexpression",
    "mtrR_promoter_ins2bp": "efflux_pump_overexpression",
    "mtrR_promoter_mosaic_1_present": "efflux_pump_overexpression",
    "mtrR_promoter_mosaic_2_present": "efflux_pump_overexpression",
    "mtrR_promoter_mosaic_3_present": "efflux_pump_overexpression",
    "mtrD_mosaic_1_present":         "efflux_pump_overexpression",
    "mtrD_mosaic_2_present":         "efflux_pump_overexpression",
    "mtrD_mosaic_3_present":         "efflux_pump_overexpression",
    "mtrD_mosaic_ambiguous_present": "efflux_pump_overexpression",

    "macAB_promoter_mut": "moderate_azithromycin_resistance",
    "norM_promoter_mut":  "norM_efflux_upregulation",

    "rpsJ_V57M": "tetracycline_chromosomal_resistance",

    "folP_ins166": "sulfonamide_resistance",
    "folP_ins167": "sulfonamide_resistance",
    "folP_ins168": "sulfonamide_resistance",
    "folP_ins169": "sulfonamide_resistance",
    "folP_ins170": "sulfonamide_resistance",
    "folP_ins171": "sulfonamide_resistance",
    "folP_ins172": "sulfonamide_resistance",
    "folP_F31L":  "sulfonamide_resistance",
    "folP_G190D": "sulfonamide_resistance",
    "folP_G190R": "sulfonamide_resistance",
    "folP_P192S": "sulfonamide_resistance",
    "folP_R228S": "sulfonamide_resistance",

    "rplD_K51E": "moderate_azithromycin_resistance",
    "rplD_Q66H": "moderate_azithromycin_resistance",
    "rplD_Q66K": "moderate_azithromycin_resistance",
    "rplD_G68D": "moderate_azithromycin_resistance",
    "rplD_G68C": "moderate_azithromycin_resistance",
    "rplD_G70D": "moderate_azithromycin_resistance",

    "rpsE_T24P":   "spectinomycin_resistance",
    "rpsE_delV25": "spectinomycin_resistance",
    "rpsE_delK26": "spectinomycin_resistance",
    "rpsE_delK25": "spectinomycin_resistance",
    "rpsE_delV26": "spectinomycin_resistance",
    "rpsE_delV27": "spectinomycin_resistance",
    "rpsE_K28E":   "spectinomycin_resistance",
    "rpsE_C192R":  "spectinomycin_resistance",
    "rpsE_W197L":  "spectinomycin_resistance",

    "ponA_L421P": "reduced_penicillin_susceptibility",

    "porB_G120K": "reduced_beta_lactam_susceptibility",
    "porB_A121D": "reduced_beta_lactam_susceptibility",
    "porB_A121N": "reduced_beta_lactam_susceptibility",
    "porB_G120D": "reduced_beta_lactam_susceptibility",

    "rpoB_P157L": "reduced_beta_lactam_susceptibility",
    "rpoB_G158V": "reduced_beta_lactam_susceptibility",
    "rpoB_R201H": "ceftriaxone_reduced_susceptibility",

    "rpoD_delD92": "reduced_beta_lactam_susceptibility",
    "rpoD_delD93": "reduced_beta_lactam_susceptibility",
    "rpoD_delD94": "reduced_beta_lactam_susceptibility",
    "rpoD_delA95": "reduced_beta_lactam_susceptibility",
    "rpoD_E98K":   "reduced_beta_lactam_susceptibility",

    "pilQ_E666K": "reduced_beta_lactam_susceptibility",
}

_POSITION_FILTER_GENES: frozenset[str] = frozenset({"penA", "porB", "mtrR", "gyrA"})

_PENA_MOSAIC_ASSOCIATED: frozenset[str] = frozenset({
    "penA_A311V", "penA_I312M", "penA_V316T", "penA_V316P",
    "penA_T483S", "penA_F504L", "penA_A501V", "penA_A501P", "penA_A501T",
    "penA_G542S", "penA_G545S", "penA_N512Y", "penA_T534A",
    "penA_A510V", "penA_A516G", "penA_H541N", "penA_P551S", "penA_P551L",
    "penA_P552V", "penA_K555Q", "penA_I556V", "penA_I566V", "penA_A574V",
    "penA_ins345", "penA_ins346", "penA_ins573",
})
_PENA_MOSAIC_ALLELE_THRESHOLD = 3

def _build_cdc_positions() -> dict[str, set[int]]:
    positions: dict[str, set[int]] = {}
    for key in CDC_RULES:
        gene, _, mut = key.partition("_")
        if gene not in _POSITION_FILTER_GENES:
            continue
        m = _re.search(r"(\d+)", mut)
        if m:
            positions.setdefault(gene, set()).add(int(m.group(1)))
    return positions

_CDC_GENE_POSITIONS: dict[str, set[int]] = _build_cdc_positions()
_CDC_GENES: frozenset[str] = frozenset(k.split("_")[0] for k in CDC_RULES)

_MTRR_EFFLUX_KEYS: frozenset[str] = frozenset({
    "mtrR_A39T", "mtrR_G45D",
    "mtrR_promoter_del35A", "mtrR_promoter_AtoC", "mtrR_promoter_mtr120", "mtrR_promoter_ins2bp",
    "mtrR_promoter_mosaic_1_present", "mtrR_promoter_mosaic_2_present", "mtrR_promoter_mosaic_3_present",
    "mtrD_mosaic_1_present", "mtrD_mosaic_2_present", "mtrD_mosaic_3_present",
    "mtrD_mosaic_ambiguous_present",
})

# Genes where any premature stop codon confers a known phenotype.
_TRUNCATION_RULES = {
    "mtrR": "efflux_pump_overexpression",  # loss of repressor, mtrCDE overexpression
}

AMR_TIER: dict[str, int] = {
    "penA_A501V": 1, "penA_A501P": 1, "penA_A501T": 1,
    "penA_G545S": 1, "penA_N512Y": 1,
    "penA_P551S": 1, "penA_P551L": 1,
    "rpoB_R201H": 1,
    "23SrRNA_A2045G": 1,
    "gyrA_A67P": 1, "gyrA_S91F": 1, "gyrA_S91Y": 1,
    "gyrA_D95G": 1, "gyrA_D95N": 1, "gyrA_D95A": 1, "gyrA_D95E": 1,
    "parC_D86N": 1, "parC_S87R": 1, "parC_S87W": 1, "parC_A87G": 1,
    "16SrRNA_C1192T": 1,
    "rpsE_T24P": 1, "rpsE_delV25": 1, "rpsE_delK26": 1,
    "rpsE_delK25": 1, "rpsE_delV26": 1, "rpsE_delV27": 1, "rpsE_K28E": 1,
    "blaTEM-1": 1, "blaTEM-135": 1,
    "tet-M": 1,
    "ermA": 1, "ermB": 1, "ermC": 1, "ermF": 1,
    "ereA": 1, "ereB": 1, "mef": 1,
    "aac-aph": 1,
    "gyrB_G447V": 3, "gyrB_D429N": 3, "gyrB_A508V": 3, "gyrB_K450T": 3, "gyrB_E466D": 3,
    "parE_G410V": 3, "parE_D420N": 3, "parE_L445H": 3, "parE_T451I": 3, "parE_T451S": 3,
    "norM_promoter_mut": 3,
}

PLASMID_CDC_RULES = {
    "blaTEM-1":   "penicillin_resistance",
    "blaTEM-135": "penicillin_resistance",
    "tet-M":    "tetracycline_resistance",
    "ermA":     "macrolide_resistance",
    "ermB":     "macrolide_resistance",
    "ermC":     "macrolide_resistance",
    "ermF":     "macrolide_resistance",
    "ereA":     "macrolide_resistance",
    "ereB":     "macrolide_resistance",
    "mef":      "macrolide_resistance",
    "aac-aph":  "aminoglycoside_resistance",
}

# Maps each CDC phenotype to the antibiotic class 
_PHENOTYPE_TO_CLASS: dict[str, str] = {
    "ceftriaxone_reduced_susceptibility":  "ceftriaxone",
    "reduced_beta_lactam_susceptibility":  "beta_lactam",
    "penicillin_resistance":               "penicillin",
    "reduced_penicillin_susceptibility":   "penicillin",
    "high_level_azithromycin_resistance":  "azithromycin",
    "moderate_azithromycin_resistance":    "azithromycin",
    "macrolide_resistance":                "azithromycin",
    "ciprofloxacin_resistant":             "fluoroquinolone",
    "ciprofloxacin_intermediate_parC":     "fluoroquinolone",
    "fluoroquinolone_minor_parE":          "fluoroquinolone",
    "fluoroquinolone_minor_gyrB":          "fluoroquinolone",
    "tetracycline_resistance":             "tetracycline",
    "tetracycline_chromosomal_resistance": "tetracycline",
    "sulfonamide_resistance":              "sulfonamide",
    "spectinomycin_resistance":            "spectinomycin",
    "aminoglycoside_resistance":           "aminoglycoside",
    "efflux_pump_overexpression":          "beta_lactam",
    "efflux_tetracycline_contribution":    "tetracycline",
    "macrolide_efflux_upregulation":       "azithromycin",
    "norM_efflux_upregulation":            "fluoroquinolone",
}

_HIGH_IMPACT = frozenset({
    "ceftriaxone_reduced_susceptibility",
    "high_level_azithromycin_resistance",
})

_MODERATE_IMPACT = frozenset({
    "ciprofloxacin_resistant",
    "ciprofloxacin_intermediate_parC",
    "penicillin_resistance",
    "tetracycline_resistance",
    "tetracycline_chromosomal_resistance",
    "sulfonamide_resistance",
    "spectinomycin_resistance",
    "macrolide_resistance",
    "moderate_azithromycin_resistance",
    "aminoglycoside_resistance",
    "reduced_penicillin_susceptibility",
})

_LOW_IMPACT = frozenset({
    "efflux_pump_overexpression",
    "efflux_tetracycline_contribution",
    "reduced_beta_lactam_susceptibility",
    "fluoroquinolone_minor_gyrB",
    "fluoroquinolone_minor_parE",
    "norM_efflux_upregulation",
})

@lru_cache(maxsize=None)
def _load_freq_panel_raw() -> dict:
    if not _FREQ_PANEL_PATH.exists():
        return {}
    try:
        with open(_FREQ_PANEL_PATH) as f:
            return json.load(f)
    except Exception:
        return {}


def load_synonymous_panel() -> dict:
    return _load_freq_panel_raw().get("frequencies", {})


def synonymous_panel_meta() -> dict:
    return _load_freq_panel_raw().get("meta", {"n_genomes": 0, "genomes": []})


def _significant_classes(phenotypes: list) -> set:
    return {
        _PHENOTYPE_TO_CLASS[p] for p in phenotypes
        if p in _PHENOTYPE_TO_CLASS and p not in _LOW_IMPACT
    }


def count_resistance_classes(phenotypes: list) -> int:
    """Return the number of distinct antibiotic classes with resistance."""
    return len(_significant_classes(phenotypes))


def classify_resistance_category(phenotypes: list) -> str:
    """
    Assign a clinical resistance category from the detected phenotypes.

    Priority, highest is first:
      MDR              — ≥2 antibiotic classes affected
      high_resistance  — ≥1 high-impact phenotype (ceftriaxone, HLAzithro)
      moderate_resistance — ≥1 moderate-impact phenotype
      low_resistance   — only minor mechanisms (efflux, reduced susceptibility)
      susceptible      — no resistance
    """
    pheno_set = set(phenotypes)
    if not pheno_set or pheno_set == {"wildtype"}:
        return "susceptible"

    n_classes = count_resistance_classes(phenotypes)
    if n_classes >= 2:
        return "MDR"
    if pheno_set & _HIGH_IMPACT:
        return "high_resistance"
    if pheno_set & _MODERATE_IMPACT:
        return "moderate_resistance"
    if pheno_set & _LOW_IMPACT:
        return "low_resistance"
    return "susceptible"


# WHO reference strain profiles for MDR/XDR identification
# Source: Unemo et al. 2016 (J Antimicrob Chemother 71:3096-3108); WHO 2012 and 2016 panels.

WHO_STRAIN_PROFILES: dict[str, dict] = {
    "WHO F":  {"mdr_class": "susceptible", "resistance_classes": set()},
    "WHO G":  {"mdr_class": "MDR", "resistance_classes": {"fluoroquinolone", "tetracycline", "ceftriaxone"}},
    "WHO K":  {"mdr_class": "MDR", "resistance_classes": {"azithromycin", "fluoroquinolone", "tetracycline"}},
    "WHO L":  {"mdr_class": "MDR", "resistance_classes": {"fluoroquinolone", "tetracycline", "penicillin", "ceftriaxone"}},
    "WHO M":  {"mdr_class": "MDR", "resistance_classes": {"fluoroquinolone", "tetracycline", "ceftriaxone"}},
    "WHO N":  {"mdr_class": "MDR", "resistance_classes": {"fluoroquinolone", "ceftriaxone"}},
    "WHO O":  {"mdr_class": "MDR", "resistance_classes": {"fluoroquinolone", "penicillin"}},
    "WHO P":  {"mdr_class": "XDR", "resistance_classes": {"ceftriaxone", "azithromycin", "fluoroquinolone", "tetracycline", "penicillin", "sulfonamide"}},
    "WHO U":  {"mdr_class": "MDR", "resistance_classes": {"fluoroquinolone", "tetracycline", "ceftriaxone"}},
    "WHO V":  {"mdr_class": "MDR", "resistance_classes": {"fluoroquinolone", "tetracycline"}},
    "WHO W":  {"mdr_class": "MDR", "resistance_classes": {"fluoroquinolone", "tetracycline"}},
    "WHO X":  {"mdr_class": "MDR", "resistance_classes": {"ceftriaxone"}},
    "WHO Y":  {"mdr_class": "MDR", "resistance_classes": {"azithromycin", "fluoroquinolone", "tetracycline"}},
    "WHO Z":  {"mdr_class": "MDR", "resistance_classes": {"ceftriaxone", "fluoroquinolone"}},
}


def match_who_strain(phenotypes: list) -> list[dict]:
    """
    Compare detected phenotypes against WHO reference strain profiles.
    A match is reported when the sample's resistance classes fully cover the
    WHO strain's defining classes (sample may be equally or more resistant).
    Returns matches sorted by severity: XDR first, then MDR.
    """
    detected_classes = {
        _PHENOTYPE_TO_CLASS[p] for p in phenotypes if p in _PHENOTYPE_TO_CLASS
    }
    matches = []
    for strain, profile in WHO_STRAIN_PROFILES.items():
        required = profile["resistance_classes"]
        if not required:
            continue
        if required.issubset(detected_classes):
            matches.append({
                "strain":          strain,
                "mdr_class":       profile["mdr_class"],
                "matched_classes": sorted(required),
            })
    order = {"XDR": 0, "MDR": 1}
    matches.sort(key=lambda x: order.get(x["mdr_class"], 2))
    return matches


# Clinical impact tier per mutation, derived from the same CDC_RULES phenotype
# mapping and _HIGH_IMPACT/_MODERATE_IMPACT/_LOW_IMPACT classification used for
# resistance-category scoring (see estimate_failure_probability). No standalone
# numeric weights: replaces the earlier per-mutation FAILURE_WEIGHTS table, whose
# individual decimal values (e.g. penA_A501P=0.55 vs penA_G545S=0.20) had no
# literature source beyond the coarse tier a mutation already belongs to.

def mutation_impact(gene: str, mutation: str | None = None) -> str | None:
    """Return 'high' / 'moderate' / 'low' clinical impact tier for a mutation
    (or a plasmid gene, when mutation is None), or None if unclassified."""
    phenotype = PLASMID_CDC_RULES.get(gene) if mutation is None else CDC_RULES.get(f"{gene}_{mutation}")
    if phenotype in _HIGH_IMPACT:
        return "high"
    if phenotype in _MODERATE_IMPACT:
        return "moderate"
    if phenotype in _LOW_IMPACT:
        return "low"
    return None

_STD_ALTERNATIVES = [
    {"condition": "beta-lactam allergy", "regimen": "spectinomycin 2g IM + azithromycin 2g orally (single dose)"},
    {"condition": "beta-lactam allergy (alternative)", "regimen": "gentamicin 240mg IM + azithromycin 2g orally (single dose)"},
    {"condition": "IM route contraindicated", "regimen": "cefixime 400mg + azithromycin 2g orally (single dose)"},
    {"condition": "fluoroquinolone susceptibility confirmed", "regimen": "ciprofloxacin 500mg orally (single dose) — avoid in pregnancy and age >60y"},
]

THERAPY_RULES = {
    "ceftriaxone_reduced_susceptibility": {
        "avoid": ["ceftriaxone", "cefixime"],
        "recommend": ["gentamicin 240mg IM + azithromycin 2g orally (single dose)"],
        "alternatives": [
            {"condition": "gentamicin unavailable", "regimen": "spectinomycin 2g IM + azithromycin 2g orally (single dose)"},
            {"condition": "specialist management", "regimen": "three-site culture + AMR testing required; notify public health authorities"},
        ],
    },
    "high_level_azithromycin_resistance": {
        "avoid": ["azithromycin"],
        "recommend": ["ceftriaxone 1g IM (single dose)"],
        "alternatives": [
            {"condition": "beta-lactam allergy", "regimen": "spectinomycin 2g IM (single dose) — specialist review required; no azithromycin"},
            {"condition": "beta-lactam allergy (alternative)", "regimen": "gentamicin 240mg IM (single dose) — specialist review required"},
            {"condition": "IM route contraindicated", "regimen": "cefixime 400mg orally (single dose) — TOC mandatory"},
            {"condition": "fluoroquinolone susceptibility confirmed", "regimen": "ciprofloxacin 500mg orally (single dose) — avoid in pregnancy"},
        ],
    },
    "moderate_azithromycin_resistance": {
        "avoid": [],
        "recommend": ["ceftriaxone 1g IM + azithromycin 2g orally (single dose — monitor MIC)"],
        "alternatives": _STD_ALTERNATIVES,
    },
    "ciprofloxacin_resistant": {
        "avoid": ["ciprofloxacin"],
        "recommend": ["ceftriaxone 1g IM + azithromycin 2g orally (single dose)"],
        "alternatives": [
            {"condition": "beta-lactam allergy", "regimen": "spectinomycin 2g IM + azithromycin 2g orally (single dose)"},
            {"condition": "beta-lactam allergy (alternative)", "regimen": "gentamicin 240mg IM + azithromycin 2g orally (single dose)"},
            {"condition": "IM route contraindicated", "regimen": "cefixime 400mg + azithromycin 2g orally (single dose)"},
        ],
    },
    "ciprofloxacin_intermediate_parC": {
        "avoid": ["ciprofloxacin"],
        "recommend": ["ceftriaxone 1g IM + azithromycin 2g orally (single dose)"],
        "alternatives": [
            {"condition": "beta-lactam allergy", "regimen": "spectinomycin 2g IM + azithromycin 2g orally (single dose)"},
            {"condition": "beta-lactam allergy (alternative)", "regimen": "gentamicin 240mg IM + azithromycin 2g orally (single dose)"},
            {"condition": "IM route contraindicated", "regimen": "cefixime 400mg + azithromycin 2g orally (single dose)"},
        ],
    },
    "penicillin_resistance": {
        "avoid": ["penicillin", "ampicillin"],
        "recommend": ["ceftriaxone 1g IM + azithromycin 2g orally (single dose)"],
        "alternatives": _STD_ALTERNATIVES,
    },
    "tetracycline_resistance": {
        "avoid": ["tetracycline", "doxycycline"],
        "recommend": ["ceftriaxone 1g IM + azithromycin 2g orally (single dose)"],
        "alternatives": _STD_ALTERNATIVES,
    },
    "sulfonamide_resistance": {
        "avoid": ["sulfonamides"],
        "recommend": ["ceftriaxone 1g IM + azithromycin 2g orally (single dose)"],
        "alternatives": _STD_ALTERNATIVES,
    },
    "spectinomycin_resistance": {
        "avoid": ["spectinomycin"],
        "recommend": ["ceftriaxone 1g IM + azithromycin 2g orally (single dose)"],
        "alternatives": [
            {"condition": "beta-lactam allergy", "regimen": "gentamicin 240mg IM + azithromycin 2g orally (single dose)"},
            {"condition": "IM route contraindicated", "regimen": "cefixime 400mg + azithromycin 2g orally (single dose)"},
            {"condition": "fluoroquinolone susceptibility confirmed", "regimen": "ciprofloxacin 500mg orally (single dose) — avoid in pregnancy"},
        ],
    },
    "macrolide_resistance": {
        "avoid": ["azithromycin"],
        "recommend": ["ceftriaxone 1g IM (single dose)"],
        "alternatives": [
            {"condition": "beta-lactam allergy", "regimen": "spectinomycin 2g IM (single dose) — specialist review required; no azithromycin"},
            {"condition": "beta-lactam allergy (alternative)", "regimen": "gentamicin 240mg IM (single dose) — specialist review required"},
            {"condition": "IM route contraindicated", "regimen": "cefixime 400mg orally (single dose) — TOC mandatory"},
            {"condition": "fluoroquinolone susceptibility confirmed", "regimen": "ciprofloxacin 500mg orally (single dose) — avoid in pregnancy"},
        ],
    },
    "aminoglycoside_resistance": {
        "avoid": ["gentamicin"],
        "recommend": ["ceftriaxone 1g IM + azithromycin 2g orally (single dose)"],
        "alternatives": _STD_ALTERNATIVES,
    },
    "reduced_beta_lactam_susceptibility": {
        "avoid": [],
        "recommend": ["ceftriaxone 1g IM + azithromycin 2g orally (single dose — monitor MIC)"],
        "alternatives": _STD_ALTERNATIVES,
    },
    "fluoroquinolone_minor_gyrB": {
        "avoid": [],
        "recommend": ["ceftriaxone 1g IM + azithromycin 2g orally (single dose — gyrB minor contributor; significant only with gyrA mutations)"],
        "alternatives": _STD_ALTERNATIVES,
    },
    "fluoroquinolone_minor_parE": {
        "avoid": [],
        "recommend": ["ceftriaxone 1g IM + azithromycin 2g orally (single dose — parE minor contributor; significant only with gyrA/parC mutations)"],
        "alternatives": _STD_ALTERNATIVES,
    },
    "efflux_pump_overexpression": {
        "avoid": [],
        "recommend": ["ceftriaxone 1g IM + azithromycin 2g orally (single dose — monitor MIC)"],
        "alternatives": _STD_ALTERNATIVES,
    },
    "tetracycline_chromosomal_resistance": {
        "avoid": ["tetracycline", "doxycycline"],
        "recommend": ["ceftriaxone 1g IM + azithromycin 2g orally (single dose)"],
        "alternatives": _STD_ALTERNATIVES,
    },
    "reduced_penicillin_susceptibility": {
        "avoid": ["penicillin", "ampicillin"],
        "recommend": ["ceftriaxone 1g IM + azithromycin 2g orally (single dose)"],
        "alternatives": _STD_ALTERNATIVES,
    },
    "wildtype": {
        "avoid": [],
        "recommend": ["ceftriaxone 1g IM + azithromycin 2g orally (single dose)"],
        "alternatives": _STD_ALTERNATIVES,
    },
}


def run_cmd(cmd, input_data=None):
    proc = subprocess.run(
        cmd,
        input=input_data,
        capture_output=True,
        text=True
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\n"
            f"STDERR:\n{proc.stderr}"
        )
    return proc.stdout


def _run_cmd_binary(cmd: list) -> bytes:
    proc = subprocess.run(cmd, capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\nSTDERR:\n{proc.stderr.decode()}"
        )
    return proc.stdout


def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)

# AMR functions

def detect_presence(ref_gene: Path, contigs: Path,
                    min_pid: float = 90.0, min_cov: float = 0.80) -> bool:
    proc = subprocess.run(
        [
            BLASTN,
            "-query",           str(ref_gene),
            "-subject",         str(contigs),
            "-outfmt",          "6 pident length qlen",
            "-perc_identity",   str(min_pid),
            "-max_target_seqs", "1",
            "-dust",            "no",
        ],
        capture_output=True, text=True,
    )
    for line in proc.stdout.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        try:
            cov = int(parts[1]) / int(parts[2])
        except (ValueError, ZeroDivisionError):
            continue
        if cov >= min_cov:
            return True
    return False


def align_and_sort(ref_gene: Path, contigs: Path, outdir: Path):
    sam = outdir / f"{ref_gene.stem}.sam"
    bam = outdir / f"{ref_gene.stem}.bam"
    sorted_bam = outdir / f"{ref_gene.stem}.sorted.bam"

    # asm20: assembled contigs vs reference, tolerates up to ~20% divergence
    sam_out = run_cmd([MINIMAP2, "-a", "-x", "asm20", str(ref_gene), str(contigs)])
    sam.write_text(sam_out)

    bam.write_bytes(_run_cmd_binary(["samtools", "view", "-bS", str(sam)]))

    if _samtools_sort_major() >= 1:
        run_cmd(["samtools", "sort", "-o", str(sorted_bam), str(bam)])
    else:
        run_cmd(["samtools", "sort", "-f", str(bam), str(sorted_bam)])
    run_cmd(["samtools", "index", str(sorted_bam)])

    sam.unlink(missing_ok=True)
    bam.unlink(missing_ok=True)

    return sorted_bam


def call_variants(ref_gene: Path, contigs: Path, outdir: Path) -> Path:
    ensure_dir(outdir)

    fai = Path(str(ref_gene) + ".fai")
    if not fai.exists():
        with _FAIDX_LOCK:
            if not fai.exists():
                subprocess.run(["samtools", "faidx", str(ref_gene)], capture_output=True)

    sorted_bam = align_and_sort(ref_gene, contigs, outdir)
    vcf = outdir / f"{ref_gene.stem}.vcf"

    mpileup = subprocess.run(
        ["bcftools", "mpileup", "-f", str(ref_gene), str(sorted_bam)],
        capture_output=True
    )
    if mpileup.returncode != 0:
        raise RuntimeError(
            f"bcftools mpileup failed for {ref_gene.stem}\nSTDERR: {mpileup.stderr.decode()}"
        )

    call_out = subprocess.run(
        ["bcftools", "call", "-mv", "--ploidy", "1", "-o", str(vcf)],
        input=mpileup.stdout,
        capture_output=True
    )
    if call_out.returncode != 0:
        raise RuntimeError(
            f"bcftools call failed for {ref_gene.stem}\nSTDERR: {call_out.stderr.decode()}"
        )

    return vcf


def parse_vcf(vcf_path: Path, ref_gene: Path = None) -> tuple[list, list]:
    """
    Parse a VCF and return (resistance_mutations, synonymous_mutations).

    SNPs that share a codon are combined before translation so that compound
    multi-nucleotide variants (e.g. GGT→AAG = G120K from three individual SNPs)
    are reported as the correct amino-acid change rather than incorrect
    single-substitution intermediates.
    """
    is_rna = ref_gene is not None and ref_gene.stem in _RNA_GENES
    ref_seq = None
    aa_offset = _AA_POSITION_OFFSET.get(ref_gene.stem, 0) if ref_gene else 0
    if ref_gene is not None and not is_rna:
        try:
            ref_seq = _read_fasta_seq(ref_gene)
        except Exception:
            pass

    mutations = []
    synonymous = []

    gene_name = ref_gene.stem if ref_gene else None
    known_positions = _CDC_GENE_POSITIONS.get(gene_name) if gene_name in _POSITION_FILTER_GENES else None

    def _at_known_pos(mut_str: str) -> bool:
        if known_positions is None:
            return True
        m = _re.search(r"(\d+)", mut_str)
        return bool(m and int(m.group(1)) in known_positions)

    snp_rows: list[tuple[int, str, str]] = []   # (1-based pos, ref_nuc, alt_nuc)
    indel_rows: list[tuple[int, str, str]] = []

    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 5:
                continue
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            ref_nuc = fields[3]
            alt_nuc = fields[4]
            is_snp = len(ref_nuc) == 1 and len(alt_nuc) == 1
            if is_snp:
                snp_rows.append((pos, ref_nuc, alt_nuc))
            else:
                indel_rows.append((pos, ref_nuc, alt_nuc))

    for pos, ref_nuc, alt_nuc in indel_rows:
        is_insertion = len(alt_nuc) > len(ref_nuc) and len(ref_nuc) == 1
        is_deletion  = len(ref_nuc) > len(alt_nuc) and len(alt_nuc) == 1

        if is_insertion and not is_rna:
            ins_len = len(alt_nuc) - len(ref_nuc)
            if ins_len % 3 == 0:
                codon_num = (pos - 1) // 3 + 1
                mut_str = f"ins{codon_num}"
                if _at_known_pos(mut_str):
                    mutations.append(mut_str)
            continue

        if is_deletion and not is_rna:
            del_len = len(ref_nuc) - len(alt_nuc)
            if del_len % 3 == 0 and ref_seq is not None:
                first_del_0   = pos
                first_codon_0 = first_del_0 // 3
                n_codons      = del_len // 3
                ref_aas = []
                for i in range(n_codons):
                    cs    = (first_codon_0 + i) * 3
                    codon = ref_seq[cs : cs + 3]
                    if len(codon) == 3:
                        ref_aas.append(_translate_codon(codon))
                if len(ref_aas) == 1:
                    mut_str = f"del{ref_aas[0]}{first_codon_0 + 1}"
                    if _at_known_pos(mut_str):
                        mutations.append(mut_str)
                elif len(ref_aas) > 1:
                    last = first_codon_0 + n_codons
                    mut_str = f"del{ref_aas[0]}{first_codon_0 + 1}_{ref_aas[-1]}{last}"
                    if _at_known_pos(mut_str):
                        mutations.append(mut_str)

    if is_rna:
        _rna_allowed = _RNA_GENE_POSITIONS.get(ref_gene.stem) if ref_gene else None
        for pos, ref_nuc, alt_nuc in snp_rows:
            if _rna_allowed is None or pos in _rna_allowed:
                mutations.append(f"{ref_nuc}{pos}{alt_nuc}")
    elif ref_seq is not None:
        from collections import defaultdict as _dd
        by_codon: dict[int, list[tuple[int, str, str]]] = _dd(list)
        for pos, ref_nuc, alt_nuc in snp_rows:
            idx = pos - 1
            if idx < 0 or idx >= len(ref_seq):
                continue
            if ref_seq[idx].upper() != ref_nuc.upper():
                continue
            by_codon[(idx) // 3].append((pos, ref_nuc, alt_nuc))

        for codon_idx, variants in by_codon.items():
            codon_start = codon_idx * 3
            ref_codon = ref_seq[codon_start:codon_start + 3]
            if len(ref_codon) < 3:
                continue

            alt_codon_list = list(ref_codon)
            for pos, ref_nuc, alt_nuc in variants:
                codon_pos = (pos - 1) % 3
                alt_codon_list[codon_pos] = alt_nuc.upper()
            alt_codon = "".join(alt_codon_list)

            ref_aa = _translate_codon(ref_codon)
            alt_aa = _translate_codon(alt_codon)

            if ref_aa == alt_aa:
                for pos, ref_nuc, alt_nuc in variants:
                    if ref_aa not in ("*", "X"):
                        synonymous.append({
                            "nuc": f"{ref_nuc.upper()}{pos}{alt_nuc.upper()}",
                            "aa":  ref_aa,
                            "pos": codon_idx + 1 + aa_offset,
                        })
            else:
                mut_str = f"{ref_aa}{codon_idx + 1 + aa_offset}{alt_aa}"
                if _at_known_pos(mut_str):
                    mutations.append(mut_str)
    else:
        for pos, ref_nuc, alt_nuc in snp_rows:
            mutations.append(f"{ref_nuc}{pos}{alt_nuc}")

    return mutations, synonymous


def infer_cdc_phenotype(chrom_mutations: dict, plasmid_results: dict = None):
    phenotypes = []
    mosaic_hits: list[str] = []

    for gene, muts in chrom_mutations.items():
        for m in (muts or []):
            key = f"{gene}_{m}"
            if key in CDC_RULES:
                phenotypes.append(CDC_RULES[key])
                if key in _MTRR_EFFLUX_KEYS:
                    phenotypes.append("efflux_tetracycline_contribution")
                if key in _PENA_MOSAIC_ASSOCIATED:
                    mosaic_hits.append(key)
            elif m.endswith("*") and gene in _TRUNCATION_RULES:
                phenotypes.append(_TRUNCATION_RULES[gene])
                phenotypes.append("efflux_tetracycline_contribution")

    for gene, status in (plasmid_results or {}).items():
        if status == "present" and gene in PLASMID_CDC_RULES:
            phenotypes.append(PLASMID_CDC_RULES[gene])

    if len(mosaic_hits) >= _PENA_MOSAIC_ALLELE_THRESHOLD:
        phenotypes.append("penA_allele_likely_mosaic")

    return list(dict.fromkeys(phenotypes)) or ["wildtype"]


_LASTLINE_CLASSES = frozenset({"ceftriaxone", "azithromycin"})

_RESISTANCE_SCORE_BY_CATEGORY = {
    "susceptible":         0.0,
    "low_resistance":      0.2,
    "moderate_resistance": 0.4,
    "high_resistance":     0.6,
    "MDR":                 0.8,
}
_XDR_SCORE = 1.0


def _is_xdr_pattern(phenotypes: list) -> bool:
    classes = _significant_classes(phenotypes)
    if not _LASTLINE_CLASSES.issubset(classes):
        return False
    return len(classes - _LASTLINE_CLASSES) >= 2


def estimate_failure_probability(chrom_mutations: dict, plasmid_results: dict = None):
    phenotypes = infer_cdc_phenotype(chrom_mutations, plasmid_results)
    if _is_xdr_pattern(phenotypes):
        return _XDR_SCORE
    category = classify_resistance_category(phenotypes)
    return _RESISTANCE_SCORE_BY_CATEGORY.get(category, 0.0)


def recommend_therapy(phenotypes: list):
    priority = [
        "ceftriaxone_reduced_susceptibility",
        "high_level_azithromycin_resistance",
        "ciprofloxacin_resistant",
        "ciprofloxacin_intermediate_parC",
        "penicillin_resistance",
        "reduced_penicillin_susceptibility",
        "macrolide_resistance",
        "aminoglycoside_resistance",
        "tetracycline_resistance",
        "tetracycline_chromosomal_resistance",
        "sulfonamide_resistance",
        "spectinomycin_resistance",
        "moderate_azithromycin_resistance",
        "efflux_pump_overexpression",
    ]
    for p in priority:
        if p in phenotypes:
            return THERAPY_RULES[p]
    return THERAPY_RULES["wildtype"]


def detect_mosaic_pena(contigs: Path) -> dict:
    """
    Detect potential mosaic penA alleles using BLASTN alignment identity.
    Identity thresholds:
      ≥ 99% → wild-type penA class
      97–99% → variant penA (non-mosaic point mutations)
      93–97% → likely mosaic penA (allele 60-type)
      < 93% → non-gonococcal penA detected
    """
    ref = GENE_DB_CROM / "penA.fasta"
    if not ref.exists():
        return {"identity": None, "mosaic_suspected": False, "allele_class": "reference not found"}
    try:
        proc = subprocess.run(
            [
                BLASTN,
                "-query", str(ref),
                "-subject", str(contigs),
                "-outfmt", "6 pident length qlen",
                "-max_hsps", "1",
                "-perc_identity", "80",
            ],
            capture_output=True,
            text=True,
            timeout=60,
        )
        best_identity = None
        best_cov = 0.0
        for line in proc.stdout.strip().splitlines():
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            try:
                pident = float(parts[0])
                length = int(parts[1])
                qlen   = int(parts[2])
            except ValueError:
                continue
            cov = length / qlen if qlen > 0 else 0.0
            if cov > best_cov:
                best_cov     = cov
                best_identity = pident

        if best_identity is None or best_cov < 0.7:
            return {"identity": None, "mosaic_suspected": False, "allele_class": "no hit (penA absent or too divergent)"}

        if best_identity >= 99.0:
            allele_class = "Wild-type penA"
        elif best_identity >= 97.0:
            allele_class = "Variant penA (non-mosaic)"
        elif best_identity >= 93.0:
            allele_class = "Likely mosaic penA (allele 60-type)"
        else:
            allele_class = "Non-gonococcal penA detected"

        mosaic = best_identity < 97.0 and best_cov >= 0.7

        return {
            "identity":         round(best_identity, 2),
            "coverage_pct":     round(best_cov * 100, 1),
            "mosaic_suspected": mosaic,
            "allele_class":     allele_class,
        }
    except FileNotFoundError:
        return {"identity": None, "mosaic_suspected": False, "allele_class": "blastn not found"}
    except Exception as e:
        return {"identity": None, "mosaic_suspected": False, "allele_class": "error", "error": str(e)}


def detect_mtrr_promoter_deletion(contigs: Path) -> dict:
    PROBE_35_WT   = "CGATTGCACGGATAAAAAGTCTTTTTT"
    PROBE_35_DEL  = "CGATTGCACGGATAAAAGTCTTTTTT"
    PROBE_35_ATOC = "CGATTGCACGGATACAAAGTCTTTTTT"
    PROBE_35_INS2 = "CGATTGCACGGATAAAAAAAGTCTTTTTT"
    PROBE_120_WT  = "TGCAAAGCAGGTTATACCTGT"
    PROBE_120_MUT = "TGCAAAGCAGATTATACCTGT"

    try:
        contigs_seq = _read_fasta_seq(contigs)

        del35A = None
        if _find(PROBE_35_WT, contigs_seq):
            del35A = False
        elif _find(PROBE_35_DEL, contigs_seq):
            del35A = True

        atoc = None
        if _find(PROBE_35_WT, contigs_seq):
            atoc = False
        elif _find(PROBE_35_ATOC, contigs_seq):
            atoc = True

        ins2bp = None
        if _find(PROBE_35_WT, contigs_seq):
            ins2bp = False
        elif _find(PROBE_35_INS2, contigs_seq):
            ins2bp = True

        mtr120 = None
        if _find(PROBE_120_WT, contigs_seq):
            mtr120 = False
        elif _find(PROBE_120_MUT, contigs_seq):
            mtr120 = True

        mutations = []
        if del35A:
            mutations.append("del35A")
        if atoc:
            mutations.append("AtoC")
        if ins2bp:
            mutations.append("ins2bp")
        if mtr120:
            mutations.append("mtr120")

        status = "wildtype" if not mutations else "+".join(mutations)
        if del35A is None and atoc is None and ins2bp is None and mtr120 is None:
            status = "not_found"

        return {
            "mtrr_promoter":        status,
            "mtrr_promoter_del35A": del35A,
            "mtrr_promoter_AtoC":   atoc,
            "mtrr_promoter_ins2bp": ins2bp,
            "mtrr_promoter_mtr120": mtr120,
        }
    except Exception as e:
        return {"mtrr_promoter": "error", "mtrr_promoter_del35A": None,
                "mtrr_promoter_AtoC": None, "mtrr_promoter_ins2bp": None,
                "mtrr_promoter_mtr120": None, "error": str(e)}


def detect_pena_ins346(contigs: Path) -> bool:
    ref_path = GENE_DB_CROM / "penA.fasta"
    if not ref_path.exists():
        return False
    try:
        ref_seq = _read_fasta_seq(ref_path)
        ins_site = 1035
        pre  = ref_seq[ins_site - 26 : ins_site]
        post = ref_seq[ins_site : ins_site + 26]
        if len(pre) < 20 or len(post) < 20:
            return False
        contigs_seq = _read_fasta_seq(contigs)
        pattern = _re.compile(
            _re.escape(pre) + r'(?:GAT|GAC)' + _re.escape(post),
            _re.IGNORECASE,
        )
        return bool(pattern.search(contigs_seq) or pattern.search(_rc(contigs_seq)))
    except Exception:
        return False




_MTRD_MOSAIC_THRESHOLDS: dict[str, tuple[float, float]] = {
    "mtrD_mosaic_1": (97.5, 0.98),
    "mtrD_mosaic_2": (95.0, 0.98),
    "mtrD_mosaic_3": (95.0, 0.98),
}


def detect_mtr_mosaics(contigs: Path) -> dict:
    gene_dir = GENE_DB_CROM

    targets: list[tuple[str, float, float]] = [
        ("mtrR_promoter_mosaic_1", 90.0, 0.80),
        ("mtrR_promoter_mosaic_2", 90.0, 0.80),
        ("mtrR_promoter_mosaic_3", 90.0, 0.80),
    ] + [(name, pid, cov) for name, (pid, cov) in _MTRD_MOSAIC_THRESHOLDS.items()]

    def _check_one(args: tuple) -> tuple[str, bool | None]:
        name, min_pid, min_cov = args
        fasta = gene_dir / f"{name}.fasta"
        if not fasta.exists():
            return name, None
        return name, detect_presence(fasta, contigs, min_pid, min_cov)

    results: dict = {}
    with ThreadPoolExecutor(max_workers=len(targets)) as pool:
        for name, present in pool.map(_check_one, targets):
            if present is not None:
                results[name] = present
    return results


def detect_macab_norM_promoter(contigs: Path) -> dict:
    PROBE_MACA_WT  = "CAGGGACCTGCGGTAGAATCCGCTT"
    PROBE_MACA_MUT = "CAGGGACCTGCGGTATAATCCGCTT"
    PROBE_NORM_WT  = "CCCGTATCCGCCGTCTGACGGCACGG"
    PROBE_NORM_MUT = "CCCGTATCCGCCGTTTGACGGCACGG"

    try:
        contigs_seq = _read_fasta_seq(contigs)

        macab_mut = None
        if _find(PROBE_MACA_WT, contigs_seq):
            macab_mut = False
        elif _find(PROBE_MACA_MUT, contigs_seq):
            macab_mut = True

        norm_mut = None
        if _find(PROBE_NORM_WT, contigs_seq):
            norm_mut = False
        elif _find(PROBE_NORM_MUT, contigs_seq):
            norm_mut = True

        return {
            "macAB_promoter_mut": macab_mut,
            "norM_promoter_mut":  norm_mut,
        }
    except Exception as e:
        return {"macAB_promoter_mut": None, "norM_promoter_mut": None, "error": str(e)}


def _call_chromosomal_genes(contigs: Path, outdir: Path) -> tuple[dict, dict, dict]:
    ref_genes = [
        g for g in sorted(GENE_DB_CROM.glob("*.fasta"))
        if g.stem not in _SKIP_CHROM_GENES
    ]
    chromosomal: dict = {}
    synonymous: dict = {}
    debug: dict = {}
    if not ref_genes:
        return chromosomal, synonymous, debug

    def _process_one(ref_gene: Path) -> tuple[str, list, list, str | None]:
        try:
            vcf = call_variants(ref_gene, contigs, outdir)
            muts, syn_muts = parse_vcf(vcf, ref_gene=ref_gene)
            return ref_gene.stem, muts, syn_muts, None
        except Exception as e:
            return ref_gene.stem, [], [], str(e)

    n_workers = min(len(ref_genes), os.cpu_count() or 4)
    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        for stem, muts, syn_muts, err in pool.map(_process_one, ref_genes):
            chromosomal[stem] = muts
            if stem in _CDC_GENES:
                synonymous[stem] = syn_muts
            if err is not None:
                debug[f"error_{stem}"] = err

    return chromosomal, synonymous, debug


def _call_essential_genes(contigs: Path, outdir: Path) -> tuple[dict, dict, dict]:
    ref_genes = sorted(GENE_DB_ESSENTIAL.glob("*.fasta"))
    nonsynonymous: dict = {}
    synonymous: dict = {}
    debug: dict = {}
    if not ref_genes:
        return nonsynonymous, synonymous, debug

    def _process_one(ref_gene: Path) -> tuple[str, list, list, str | None]:
        try:
            vcf = call_variants(ref_gene, contigs, outdir)
            muts, syn_muts = parse_vcf(vcf, ref_gene=ref_gene)
            return ref_gene.stem, muts, syn_muts, None
        except Exception as e:
            return ref_gene.stem, [], [], str(e)

    n_workers = min(len(ref_genes), os.cpu_count() or 4)
    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        for stem, muts, syn_muts, err in pool.map(_process_one, ref_genes):
            nonsynonymous[stem] = muts
            synonymous[stem] = syn_muts
            if err is not None:
                debug[f"error_{stem}"] = err

    return nonsynonymous, synonymous, debug


def run_essential_gene_scan(sample_id: str, contigs: Path) -> dict:
    outdir = PROJECT_ROOT / "results" / "essential_genes" / sample_id
    ensure_dir(outdir)
    nonsynonymous, synonymous, debug = _call_essential_genes(contigs, outdir)
    return {
        "nonsynonymous": nonsynonymous,
        "synonymous":    synonymous,
        "debug":         debug,
    }


def _call_plasmid_genes(target: Path) -> tuple[dict, dict]:
    ref_genes = sorted(GENE_DB_PLASM.glob("*.fasta"))
    plasmid: dict = {}
    debug: dict = {}
    if not ref_genes:
        return plasmid, debug

    def _detect_one(ref_gene: Path) -> tuple[str, str, str | None]:
        try:
            present = detect_presence(ref_gene, target)
            return ref_gene.stem, "present" if present else "absent", None
        except Exception as e:
            return ref_gene.stem, "error", str(e)

    n_workers = min(len(ref_genes), os.cpu_count() or 4)
    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        for stem, status, err in pool.map(_detect_one, ref_genes):
            plasmid[stem] = status
            if err is not None:
                debug[f"error_{stem}"] = err

    return plasmid, debug


def _apply_pena_insertions(contigs: Path, chromosomal: dict) -> dict:
    debug: dict = {}
    try:
        if detect_pena_ins346(contigs):
            pena_muts = chromosomal.setdefault("penA", [])
            if "ins345" not in pena_muts and "ins346" not in pena_muts:
                pena_muts.append("ins346")
    except Exception as e:
        debug["error_pena_ins346"] = str(e)
    return debug


def _apply_mtrr_promoter(contigs: Path, chromosomal: dict) -> dict:
    debug: dict = {}
    try:
        mtrr_prom = detect_mtrr_promoter_deletion(contigs)
        debug["mtrr_promoter"] = mtrr_prom.get("mtrr_promoter")
        prom_muts = []
        if mtrr_prom.get("mtrr_promoter_del35A") is True:
            prom_muts.append("del35A")
        if mtrr_prom.get("mtrr_promoter_AtoC") is True:
            prom_muts.append("AtoC")
        if mtrr_prom.get("mtrr_promoter_ins2bp") is True:
            prom_muts.append("ins2bp")
        if mtrr_prom.get("mtrr_promoter_mtr120") is True:
            prom_muts.append("mtr120")
        chromosomal["mtrR_promoter"] = prom_muts
    except Exception as e:
        debug["error_mtrR_promoter"] = str(e)
    return debug


def _apply_mtr_mosaics(contigs: Path, chromosomal: dict) -> dict:
    debug: dict = {}
    try:
        mtr_mosaics = detect_mtr_mosaics(contigs)
        for name, present in mtr_mosaics.items():
            chromosomal[name] = ["present"] if present else []
        _mtrd_hits = [
            k for k in ("mtrD_mosaic_1", "mtrD_mosaic_2", "mtrD_mosaic_3")
            if chromosomal.get(k) == ["present"]
        ]
        if len(_mtrd_hits) > 1:
            for k in ("mtrD_mosaic_1", "mtrD_mosaic_2", "mtrD_mosaic_3"):
                chromosomal.pop(k, None)
            chromosomal["mtrD_mosaic_ambiguous"] = ["present"]
            debug["mtrD_mosaic_codetection"] = _mtrd_hits
    except Exception as e:
        debug["error_mtr_mosaics"] = str(e)
    return debug


def _apply_efflux_promoters(contigs: Path, chromosomal: dict) -> dict:
    debug: dict = {}
    try:
        efflux_prom = detect_macab_norM_promoter(contigs)
        debug["macAB_norM_promoter"] = efflux_prom
        chromosomal["macAB_promoter"] = ["mut"] if efflux_prom.get("macAB_promoter_mut") is True else []
        chromosomal["norM_promoter"]  = ["mut"] if efflux_prom.get("norM_promoter_mut")  is True else []
    except Exception as e:
        debug["error_macAB_norM_promoter"] = str(e)
    return debug


def run_amr_variant_calling(sample_id: str, contigs: Path, plasmid_contigs: Path = None):
    outdir = PROJECT_ROOT / "results" / "amr" / sample_id
    ensure_dir(outdir)
    logger.info("AMR | %s | start", sample_id)

    _plasmid_search_target = plasmid_contigs if (plasmid_contigs and plasmid_contigs.exists()) else contigs
    _using_plasmid_assembly = plasmid_contigs is not None and plasmid_contigs.exists()

    chromosomal, synonymous, dbg_chrom = _call_chromosomal_genes(contigs, outdir)
    plasmid, dbg_plasm                 = _call_plasmid_genes(_plasmid_search_target)

    dbg_pena   = _apply_pena_insertions(contigs, chromosomal)
    dbg_mtrr   = _apply_mtrr_promoter(contigs, chromosomal)
    dbg_mosaic = _apply_mtr_mosaics(contigs, chromosomal)
    dbg_efflux = _apply_efflux_promoters(contigs, chromosomal)

    debug = {
        "gene_db_crom":          str(GENE_DB_CROM),
        "gene_db_plasm":         str(GENE_DB_PLASM),
        "crom_exists":           GENE_DB_CROM.exists(),
        "plasm_exists":          GENE_DB_PLASM.exists(),
        "crom_genes":            [f.stem for f in sorted(GENE_DB_CROM.glob("*.fasta"))],
        "plasm_genes":           [f.stem for f in sorted(GENE_DB_PLASM.glob("*.fasta"))],
        "plasmid_assembly_used": _using_plasmid_assembly,
        **dbg_chrom, **dbg_plasm, **dbg_pena, **dbg_mtrr, **dbg_mosaic, **dbg_efflux,
    }

    cdc_pheno    = infer_cdc_phenotype(chromosomal, plasmid)
    failure_prob = estimate_failure_probability(chromosomal, plasmid)

    _mlst: dict = {}
    _mosaic_pena: dict = {}
    _ngstar: dict = {}

    try:
        from backend.mlst import run_mlst
        _mlst = run_mlst(contigs, minimap2=MINIMAP2)
    except Exception as e:
        _mlst = {"error": str(e), "st": None, "alleles": {}}

    try:
        _mosaic_pena = detect_mosaic_pena(contigs)
    except Exception as e:
        _mosaic_pena = {"error": str(e), "mosaic_suspected": False}

    try:
        from backend.ngstar import run_ngstar
        _ngstar = run_ngstar(contigs, minimap2=MINIMAP2)
    except Exception as e:
        _ngstar = {"error": str(e), "ST": None, "alleles": {}}

    try:
        _essential = run_essential_gene_scan(sample_id, contigs)
    except Exception as e:
        _essential = {"nonsynonymous": {}, "synonymous": {}, "debug": {"error": str(e)}}

    _phenos = [p for p in cdc_pheno if p != "wildtype"]
    logger.info("AMR | %s | done phenotypes=%s failure_prob=%.2f",
                sample_id, _phenos or ["wildtype"], failure_prob)

    from .models import AMRResult
    return AMRResult(
        sample_id=sample_id,
        chromosomal=chromosomal,
        synonymous_variants=synonymous,
        plasmid=plasmid,
        cdc_phenotypes=cdc_pheno,
        resistance_category=classify_resistance_category(cdc_pheno),
        n_resistance_classes=count_resistance_classes(cdc_pheno),
        failure_probability=failure_prob,
        therapy=recommend_therapy(cdc_pheno),
        who_matches=match_who_strain(cdc_pheno),
        mlst=_mlst,
        mosaic_pena=_mosaic_pena,
        ngstar=_ngstar,
        essential_gene_mutations=_essential["nonsynonymous"],
        essential_gene_synonymous=_essential["synonymous"],
    )
