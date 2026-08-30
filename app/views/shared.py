from __future__ import annotations

import configparser
import os
import re
import sys
import tempfile
from pathlib import Path
from typing import Optional

import altair as alt
import pandas as pd
import streamlit as st


APP_DIR = Path(__file__).resolve().parent.parent
PROJECT_ROOT = APP_DIR.parent
sys.path.append(str(PROJECT_ROOT))

from backend.paths import PROJECTS_DIR, RESULTS_DIR
from backend.models import AMRResult

CONFIG_PATH = PROJECT_ROOT / "config.ini"
UPLOAD_DIR = APP_DIR / "uploads"
RESULTS_BASE = RESULTS_DIR


run_full_pipeline = None

try:
    from backend.full_pipeline import run_full_pipeline
except Exception:
    run_full_pipeline = None

try:
    from backend.qc import run_qc_pipeline
except Exception:
    run_qc_pipeline = None

try:
    from backend.assembly import run_assembly_pipeline, parse_contigs_stats as _parse_contigs_stats_backend
except Exception:
    run_assembly_pipeline = None
    _parse_contigs_stats_backend = None

try:
    from backend.amr import (
        run_amr_variant_calling,
        mutation_impact as _mutation_impact,
        AMR_TIER as _AMR_TIER,
        CDC_RULES as _CDC_RULES,
        PLASMID_CDC_RULES as _PLASMID_CDC_RULES,
        load_synonymous_panel,
        synonymous_panel_meta,
        _HIGH_IMPACT,
        _MODERATE_IMPACT,
        _TRUNCATION_RULES,
    )
except Exception:
    run_amr_variant_calling = None
    _mutation_impact = lambda gene, mutation=None: None
    _AMR_TIER = {}
    _CDC_RULES = {}
    _PLASMID_CDC_RULES = {}
    load_synonymous_panel = lambda: {}
    synonymous_panel_meta = lambda: {"n_genomes": 0, "genomes": []}
    _HIGH_IMPACT = frozenset()
    _MODERATE_IMPACT = frozenset()
    _TRUNCATION_RULES = {}

try:
    from backend.phylogeny.run_phylogeny import build_distance_tree
    from backend.phylogeny.viewer import load_tree_newick
except Exception:
    build_distance_tree = None
    load_tree_newick = None

try:
    import backend.db as db
except Exception:
    db = None


@st.cache_data(ttl=30)
def _load_project_cached(project_name: str) -> dict:
    if db is None:
        return {"samples": {}, "phylogeny": None}
    return db.load_project(project_name)


def _resistance_color(prob: float) -> tuple[str, str]:
    """
    Maps a resistance probability 0–1.
    White = susceptible, dark red = highly resistant.
    """
    stops = [
        (0.00, (255, 255, 255), (55,  65,  81 )),
        (0.15, (254, 226, 226), (153, 27,  27 )),
        (0.35, (252, 165, 165), (127, 29,  29 )),
        (0.55, (239, 68,  68 ), (255, 255, 255)),
        (0.75, (185, 28,  28 ), (255, 255, 255)),
        (1.00, (127, 7,   36 ), (255, 255, 255)),
    ]
    p = max(0.0, min(1.0, prob))
    for i in range(len(stops) - 1):
        p0, c0, t0 = stops[i]
        p1, c1, t1 = stops[i + 1]
        if p <= p1:
            frac = (p - p0) / (p1 - p0) if p1 > p0 else 0.0
            bg = tuple(int(c0[j] + frac * (c1[j] - c0[j])) for j in range(3))
            tc = tuple(int(t0[j] + frac * (t1[j] - t0[j])) for j in range(3))
            return f"rgb({bg[0]},{bg[1]},{bg[2]})", f"rgb({tc[0]},{tc[1]},{tc[2]})"
    return "rgb(127,7,36)", "rgb(255,255,255)"


_IMPACT_COLORS = {
    "high":     ("rgb(185,28,28)",  "rgb(255,255,255)"),
    "moderate": ("rgb(252,165,165)", "rgb(127,29,29)"),
    "low":      ("rgb(254,226,226)", "rgb(153,27,27)"),
    None:       ("rgb(255,255,255)", "rgb(55,65,81)"),
}


def _impact_color(tier: str | None) -> tuple[str, str]:
    """Fixed colors for the 3-tier clinical impact classification (high/moderate/low)."""
    return _IMPACT_COLORS.get(tier, _IMPACT_COLORS[None])


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def read_config(config_path: Path) -> configparser.ConfigParser:
    config = configparser.ConfigParser()
    config.read(config_path)
    return config


def plural(n: int, word: str, plural_form: str | None = None) -> str:
    """'sample' -> 'sample' if n == 1 else 'samples' (or plural_form if given)."""
    return word if n == 1 else (plural_form or f"{word}s")


def sanitize_name(name: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
    clean = re.sub(r"_+", "_", clean).strip("._-")
    return clean or "sample"


def safe_write_bytes(path: Path, data: bytes) -> Path:
    ensure_dir(path.parent)
    with open(path, "wb") as f:
        f.write(data)
    return path


def list_subdirs(path: Path) -> list[Path]:
    if not path.exists():
        return []
    return [p for p in path.iterdir() if p.is_dir()]


def count_items(path: Path, pattern: str = "*") -> int:
    if not path.exists():
        return 0
    return len(list(path.glob(pattern)))


def read_text_safely(path: Path) -> str:
    try:
        return path.read_text(errors="ignore")
    except Exception:
        return ""


def detect_read_type(filename: str) -> str:
    name = filename.lower()
    if re.search(r"(?:^|[_.-])(r?1)(?:[_.-]|\.|$)", name):
        return "R1"
    if re.search(r"(?:^|[_.-])(r?2)(?:[_.-]|\.|$)", name):
        return "R2"
    return "INTERLEAVED"


def group_paired_end(files) -> dict[str, dict[str, Optional[object]]]:
    samples: dict[str, dict[str, Optional[object]]] = {}
    for f in files:
        ftype = detect_read_type(f.name)
        name = re.sub(r"\s*\(\d+\)(?=\.[a-z0-9.]+$)", "", f.name.lower())
        for _ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq", ".fa.gz", ".fasta.gz", ".fa", ".fasta"):
            if name.endswith(_ext):
                name = name[: -len(_ext)]
                break
        sample_id = re.sub(r"[_.-]?r?[12]$", "", name) or name
        if sample_id not in samples:
            samples[sample_id] = {"R1": None, "R2": None}
        if ftype == "R2":
            samples[sample_id]["R2"] = f
        else:
            # run_assembly_pipeline treats r2=None as single-end.
            samples[sample_id]["R1"] = f
    return samples


def save_uploaded_file(uploaded_file, target_dir: Path = UPLOAD_DIR) -> Path:
    ensure_dir(target_dir)
    filename = sanitize_name(uploaded_file.name)
    save_path = target_dir / filename
    if save_path.exists():
        stem = save_path.stem
        suffix = "".join(save_path.suffixes)
        counter = 1
        while save_path.exists():
            save_path = target_dir / f"{stem}_{counter}{suffix}"
            counter += 1
    return safe_write_bytes(save_path, uploaded_file.getbuffer())


def save_temp_file(uploaded_file) -> Path:
    """Write an uploaded file to a OS-managed temp directory.
    The caller deletes it after use."""
    suffix = "".join(Path(uploaded_file.name).suffixes) or ".tmp"
    fd, tmp_path = tempfile.mkstemp(suffix=suffix)
    try:
        os.write(fd, uploaded_file.getbuffer())
    finally:
        os.close(fd)
    return Path(tmp_path)


def cleanup_temp_files(paths: list[Path]) -> None:
    for p in paths:
        try:
            if p.is_file():
                p.unlink()
        except Exception:
            pass


def parse_contigs_stats(path_str: str) -> dict | None:
    if not path_str or _parse_contigs_stats_backend is None:
        return None
    return _parse_contigs_stats_backend(Path(path_str)) or None


def project_dir(project_name: str) -> Path:
    return PROJECTS_DIR / sanitize_name(project_name)


def project_raw_dir(project_name: str) -> Path:
    return project_dir(project_name) / "raw"


def project_results_dir(project_name: str) -> Path:
    return project_dir(project_name) / "results"


def project_qc_dir(project_name: str) -> Path:
    return project_results_dir(project_name) / "qc"


def project_assembly_dir(project_name: str) -> Path:
    return project_results_dir(project_name) / "assembly"


def project_amr_dir(project_name: str) -> Path:
    return project_results_dir(project_name) / "amr"


def project_phylogeny_dir(project_name: str) -> Path:
    return project_results_dir(project_name) / "phylogeny"


def load_results_metrics() -> tuple[int, int, int]:
    qc_dir = RESULTS_BASE / "qc"
    asm_dir = RESULTS_BASE / "assembly"
    amr_dir = RESULTS_BASE / "amr"
    return (
        count_items(qc_dir),
        count_items(asm_dir),
        count_items(amr_dir),
    )


def collect_alerts() -> list[str]:
    qc_dir = RESULTS_BASE / "qc"
    amr_dir = RESULTS_BASE / "amr"
    alerts: list[str] = []

    for log in qc_dir.glob("*/logs.txt"):
        txt = read_text_safely(log)
        if "FAIL" in txt.upper() or "poor quality" in txt.lower():
            alerts.append(f"Poor QC detected in **{log.parent.name}**")

    for rep in qc_dir.glob("*/kraken2/kraken2_report.txt"):
        txt = read_text_safely(rep)
        if not txt:
            continue
        ng_pct = 0.0
        for line in txt.splitlines():
            parts = line.strip().split("\t")
            if len(parts) >= 5 and parts[4].strip() == "485":
                try:
                    ng_pct = float(parts[0])
                except ValueError:
                    pass
        if ng_pct < 80.0:
            alerts.append(f"Possible contamination in **{rep.parent.parent.name}** — *N. gonorrhoeae* {ng_pct:.1f}%")

    for amr_json in amr_dir.glob("*/*.json"):
        txt = read_text_safely(amr_json).lower()
        if any(gene in txt for gene in ["pena", "23s", "mtrr"]):
            alerts.append(f"AMR markers detected in **{amr_json.parent.name}**")

    return alerts


def chart_for_metrics(num_qc: int, num_asm: int, num_amr: int) -> alt.Chart:
    df = pd.DataFrame(
        {
            "Category": ["QC", "Assembly", "AMR"],
            "Samples": [num_qc, num_asm, num_amr],
        }
    )
    return (
        alt.Chart(df)
        .mark_bar()
        .encode(
            x=alt.X("Category:N", title=None),
            y=alt.Y("Samples:Q", title="Samples processed"),
            tooltip=["Category", "Samples"],
        )
        .properties(height=280)
    )


def _rebuild_pipeline_results(db_data: dict) -> dict:
    """Convert load_project() output back into full_pipeline_results session_state format."""
    samples = db_data.get("samples", {})
    results: dict = {"qc": {}, "assembly": {}, "amr": {}, "phylogeny": {}}
    for sid, data in samples.items():
        if "qc" in data:
            results["qc"][sid] = data["qc"]
        if "assembly" in data:
            results["assembly"][sid] = {
                "contigs_path": data["assembly"]["contigs_path"],
                "logs": "",
            }
        if "amr" in data:
            results["amr"][sid] = data["amr"]
    phy = db_data.get("phylogeny")
    if phy:
        cr = phy.get("cluster_report", {})
        clusters = {name: info["cluster"] for name, info in cr.items() if "cluster" in info}
        results["phylogeny"] = {
            "message":        "Tree loaded from saved project.",
            "run_id":         phy.get("run_id"),
            "tree_path":      phy.get("tree_path"),
            "clusters":       clusters,
            "cluster_report": cr,
            "matrix_path":    phy.get("matrix_path"),
            "n_samples":      phy.get("n_samples"),
        }
    return results


def init_directories() -> None:
    for path in [UPLOAD_DIR, RESULTS_BASE, PROJECTS_DIR]:
        ensure_dir(path)
    for sub in [RESULTS_BASE / "qc", RESULTS_BASE / "assembly", RESULTS_BASE / "amr"]:
        ensure_dir(sub)


def get_available_projects() -> list[str]:
    return [p.name for p in list_subdirs(PROJECTS_DIR)]


def find_trimmed_reads(project_name: str) -> list[Path]:
    qc_root = project_qc_dir(project_name)
    if not qc_root.exists():
        return []
    candidates: list[Path] = []
    for ext in ("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"):
        candidates.extend(qc_root.rglob(ext))
    return sorted(candidates, key=lambda p: ("trim" not in p.name.lower(), p.name.lower()))


def find_contigs(project_name: str) -> list[Path]:
    asm_root = project_assembly_dir(project_name)
    if not asm_root.exists():
        return []
    candidates: list[Path] = []
    for ext in ("*.fasta", "*.fa", "*.fna"):
        candidates.extend(asm_root.rglob(ext))
    return sorted(candidates, key=lambda p: p.name.lower())


_SITE_OPTIONS = ["", "urethral", "rectal", "pharyngeal", "ocular", "conjunctival", "other"]
_SEX_OPTIONS  = ["", "male", "female", "unknown"]
_META_COLS    = ["sample_id", "collection_date", "anatomical_site", "sex", "age",
                 "country", "city", "health_unit"]


_MATRIX_COLS: list[tuple] = [
    ("AZM", "23S A2045G",  "chrom",    "23SrRNA",        "A2045G"),
    ("AZM", "23S C2597T",  "chrom",    "23SrRNA",        "C2597T"),
    ("AZM", "mtrR G45D",   "chrom",    "mtrR",           "G45D"),
    ("AZM", "mtrR A39T",   "chrom",    "mtrR",           "A39T"),
    ("AZM", "mtrR -57del", "promoter", "mtrR_promoter",  "del35A"),
    ("AZM", "mtrR -56A>C", "promoter", "mtrR_promoter",  "AtoC"),
    ("AZM", "mtrR g-131a", "promoter", "mtrR_promoter",  "mtr120"),
    ("AZM", "mtrR ins2bp", "promoter", "mtrR_promoter",  "ins2bp"),
    ("AZM", "mtrR trunc",  "truncation","mtrR",           None),
    ("AZM", "mtrD mosaic", "chrom",    "mtrD_mosaic_1",  "present"),
    ("AZM", "macAB prmt",  "promoter", "macAB_promoter", "mut"),
    ("AZM", "rplD G68D",  "chrom",    "rplD",           "G68D"),
    ("AZM", "rplD G68C",  "chrom",    "rplD",           "G68C"),
    ("AZM", "rplD G70D",  "chrom",    "rplD",           "G70D"),
    ("AZM", "ermB",        "plasmid",  "ermB",           None),
    ("AZM", "ermC",        "plasmid",  "ermC",           None),
    ("AZM", "mef",         "plasmid",  "mef",            None),
    ("AZM", "ereA",        "plasmid",  "ereA",           None),
    ("AZM", "ereB",        "plasmid",  "ereB",           None),
    ("CRO", "penA I312M",  "chrom",    "penA",           "I312M"),
    ("CRO", "penA V316T",  "chrom",    "penA",           "V316T"),
    ("CRO", "penA V316P",  "chrom",    "penA",           "V316P"),
    ("CRO", "penA T483S",  "chrom",    "penA",           "T483S"),
    ("CRO", "penA A501P",  "chrom",    "penA",           "A501P"),
    ("CRO", "penA A501T",  "chrom",    "penA",           "A501T"),
    ("CRO", "penA A501V",  "chrom",    "penA",           "A501V"),
    ("CRO", "penA G542S",  "chrom",    "penA",           "G542S"),
    ("CRO", "penA G545S",  "chrom",    "penA",           "G545S"),
    ("CRO", "penA P551S",  "chrom",    "penA",           "P551S"),
    ("CRO", "penA ins345", "chrom",    "penA",           "ins345"),
    ("CRO", "penA ins346", "chrom",    "penA",           "ins346"),
    ("CRO", "penA mosaic", "mosaic",   None,             None),
    ("CRO", "rpoB P157L", "chrom",    "rpoB",           "P157L"),
    ("CRO", "rpoB G158V", "chrom",    "rpoB",           "G158V"),
    ("CRO", "rpoB R201H", "chrom",    "rpoB",           "R201H"),
    ("CRO", "rpoD E98K",  "chrom",    "rpoD",           "E98K"),
    ("CRO", "rpoD delD92","chrom",    "rpoD",           "delD92"),
    ("CRO", "rpoD delD93","chrom",    "rpoD",           "delD93"),
    ("CRO", "rpoD delD94","chrom",    "rpoD",           "delD94"),
    ("CRO", "rpoD delA95","chrom",    "rpoD",           "delA95"),
    ("PEN", "blaTEM-1",    "plasmid",  "blaTEM-1",       None),
    ("PEN", "ponA L421P",  "chrom",    "ponA",           "L421P"),
    ("PEN", "porB G120K",  "chrom",    "porB",           "G120K"),
    ("PEN", "porB A121N",  "chrom",    "porB",           "A121N"),
    ("PEN", "porB A121D",  "chrom",    "porB",           "A121D"),
    ("PEN", "pilQ E666K",  "chrom",    "pilQ",           "E666K"),
    ("CIP", "gyrA S91F",   "chrom",    "gyrA",           "S91F"),
    ("CIP", "gyrA S91Y",   "chrom",    "gyrA",           "S91Y"),
    ("CIP", "gyrA D95G",   "chrom",    "gyrA",           "D95G"),
    ("CIP", "gyrA D95N",   "chrom",    "gyrA",           "D95N"),
    ("CIP", "gyrA D95A",   "chrom",    "gyrA",           "D95A"),
    ("CIP", "gyrA D95E",   "chrom",    "gyrA",           "D95E"),
    ("CIP", "parC D86N",   "chrom",    "parC",           "D86N"),
    ("CIP", "parC S87R",   "chrom",    "parC",           "S87R"),
    ("CIP", "parC S87N",   "chrom",    "parC",           "S87N"),
    ("CIP", "parC S87I",   "chrom",    "parC",           "S87I"),
    ("CIP", "parC S88P",   "chrom",    "parC",           "S88P"),
    ("CIP", "parC E91K",   "chrom",    "parC",           "E91K"),
    ("CIP", "parE G410V",  "chrom",    "parE",           "G410V"),
    ("CIP", "norM prmt",   "promoter", "norM_promoter",  "mut"),
    ("TET", "tet-M",       "plasmid",  "tet-M",          None),
    ("TET", "rpsJ V57M",   "chrom",    "rpsJ",           "V57M"),
    ("SPE", "16S C1192T",  "chrom",    "16SrRNA",        "C1192T"),
    ("SPE", "rpsE T24P",   "chrom",    "rpsE",           "T24P"),
    ("SPE", "rpsE delV26", "chrom",    "rpsE",           "delV26"),
    ("SPE", "rpsE delV27", "chrom",    "rpsE",           "delV27"),
    ("SPE", "rpsE K28E",   "chrom",    "rpsE",           "K28E"),
    ("SUL", "folP R228S",  "chrom",    "folP",           "R228S"),
    ("SUL", "folP ins166", "chrom",    "folP",           "ins166"),
    ("AMG", "aac-aph",     "plasmid",  "aac-aph",        None),
]

_GROUP_COLORS: dict[str, str] = {
    "AZM": "#7c3aed",
    "CRO": "#dc2626",
    "PEN": "#2563eb",
    "CIP": "#16a34a",
    "TET": "#d97706",
    "SPE": "#0891b2",
    "SUL": "#be185d",
    "AMG": "#78716c",
}

_GROUP_LABELS: dict[str, str] = {
    "AZM": "Azithromycin",
    "CRO": "Ceftriaxone/Cefixime",
    "PEN": "Penicillin",
    "CIP": "Ciprofloxacin",
    "TET": "Tetracycline",
    "SPE": "Spectinomycin",
    "SUL": "Sulfonamides",
    "AMG": "Aminoglycosides",
}

_CHROM_GENE_SUPPRESS: set[str] = {
    "mtrD_mosaic_2", "mtrD_mosaic_3",
    "mtrR_promoter_mosaic_2", "mtrR_promoter_mosaic_3",
}
_CHROM_GENE_RENAME: dict[str, str] = {
    "mtrR_promoter":           "mtrR_prom",
    "mtrD_mosaic_1":           "mtrD_mosaic",
    "mtrD_mosaic_ambiguous":   "mtrD_mosaic",
    "mtrR_promoter_mosaic_1":  "mtrR_prom_mosaic",
    "macAB_promoter":          "macAB_prom",
    "norM_promoter":           "norM_prom",
}


def _mut_present(amr_result: AMRResult, lookup: str, gene_key, mut_key) -> bool:
    chrom = amr_result.chromosomal
    plasm = amr_result.plasmid
    if lookup in ("chrom", "promoter"):
        return mut_key in chrom.get(gene_key, [])
    if lookup == "plasmid":
        return plasm.get(gene_key) == "present"
    if lookup == "mosaic":
        return bool(amr_result.mosaic_pena.get("mosaic_suspected"))
    if lookup == "truncation":
        return any(m.endswith("*") for m in chrom.get(gene_key, []))
    return False


def _short_mut(display: str) -> str:
    """'mtrR G45D' → 'G45D', '23S A2059G' → 'A2059G', 'penA mosaic' → 'mosaic'"""
    parts = display.split()
    return parts[-1] if len(parts) > 1 else display


_EXPLICIT_RESISTANT_PHENOTYPES: frozenset[str] = frozenset({
    "penicillin_resistance",
    "ciprofloxacin_resistant",
    "tetracycline_resistance",
    "tetracycline_chromosomal_resistance",
    "sulfonamide_resistance",
    "spectinomycin_resistance",
    "macrolide_resistance",
    "aminoglycoside_resistance",
})


def _matrix_entry_phenotype(lookup: str, gene_key, mut_key) -> str | None:
    """The CDC phenotype label a given _MATRIX_COLS entry maps to, if any."""
    if lookup in ("chrom", "promoter"):
        return _CDC_RULES.get(f"{gene_key}_{mut_key}")
    if lookup == "plasmid":
        return _PLASMID_CDC_RULES.get(gene_key)
    if lookup == "truncation":
        return _TRUNCATION_RULES.get(gene_key)
    if lookup == "mosaic":
        return None  # identity-based detect_mosaic_pena(), independent of infer_cdc_phenotype()
    return None


def _matrix_entry_determinant(lookup: str, gene_key, mut_key) -> tuple[str, str | None]:
    """(gene, mutation) for a detected _MATRIX_COLS entry; mutation is None when it must render standalone, not slash-grouped."""
    if lookup == "plasmid":
        return (str(gene_key), None)
    if lookup == "mosaic":
        return ("penA mosaic", None)
    if lookup == "truncation":
        return (f"{gene_key}_disrupted", None)
    if lookup == "promoter":
        pw_name = {"del35A": "a-57del", "AtoC": "-56a>c", "mtr120": "g-131a"}.get(str(mut_key), str(mut_key))
        return (str(gene_key), pw_name)
    return (str(gene_key), str(mut_key))


def _join_determinants(pairs: list[tuple[str, str | None]]) -> str:
    """Group same-gene point mutations with '/', everything else with '; ' (Pathogenwatch style)."""
    groups: dict[str, list[str]] = {}
    order: list[str] = []
    standalone: list[str] = []
    for gene, mut in pairs:
        if mut is None:
            standalone.append(gene)
            continue
        if gene not in groups:
            groups[gene] = []
            order.append(gene)
        groups[gene].append(mut)
    parts = [f"{gene}_" + "/".join(groups[gene]) for gene in order] + standalone
    return "; ".join(parts)


def agent_resistance_summary(amr_result: AMRResult) -> list[dict]:
    """Per-agent Agent / Inferred resistance / Known determinants table, Pathogenwatch-style."""
    groups = [g for g in dict.fromkeys(disp[0] for disp in _MATRIX_COLS)]
    rows = []
    for group in groups:
        entries = [e for e in _MATRIX_COLS if e[0] == group]
        determinant_pairs: list[tuple[str, str | None]] = []
        phenos_found: set[str] = set()
        for _, _display, lookup, gene_key, mut_key in entries:
            if _mut_present(amr_result, lookup, gene_key, mut_key):
                determinant_pairs.append(_matrix_entry_determinant(lookup, gene_key, mut_key))
                p = _matrix_entry_phenotype(lookup, gene_key, mut_key)
                if p:
                    phenos_found.add(p)

        if phenos_found & (_HIGH_IMPACT | _EXPLICIT_RESISTANT_PHENOTYPES):
            call = "Resistant"
        elif phenos_found & _MODERATE_IMPACT:
            call = "Intermediate"
        else:
            call = "None"

        rows.append({
            "Agent":                _GROUP_LABELS.get(group, group),
            "Inferred resistance":  call,
            "Known determinants":   _join_determinants(determinant_pairs) if determinant_pairs else "-",
        })
    return rows


def render_agent_resistance_table(amr_result: AMRResult) -> None:
    """Render the Pathogenwatch-style per-agent summary table for one sample."""
    rows = agent_resistance_summary(amr_result)
    df = pd.DataFrame(rows)

    def _cell_style(v):
        if v == "Resistant":
            return "background-color:#fef2f2;color:#dc2626;font-weight:700;"
        if v == "Intermediate":
            return "background-color:#fffbeb;color:#b45309;font-weight:700;"
        return "color:#9ca3af;"

    st.dataframe(
        df.style.applymap(_cell_style, subset=["Inferred resistance"]),
        use_container_width=True,
        hide_index=True,
    )


def render_amr_mutation_matrix(results_dict: dict) -> None:
    from itertools import groupby as _groupby
    valid = {sid: r for sid, r in results_dict.items() if isinstance(r, AMRResult)}
    if not valid:
        return

    _col_order  = {disp: i for i, (_, disp, *_) in enumerate(_MATRIX_COLS)}
    _col_group  = {disp: g for g, disp, *_ in _MATRIX_COLS}
    _col_lookup = {disp: lk for _, disp, lk, *_ in _MATRIX_COLS}

    def _is_plasmid(lk: str) -> bool:
        return lk == "plasmid"

    # Presence matrix
    presence: dict[str, dict[str, bool]] = {}
    for sid, res in valid.items():
        presence[sid] = {}
        for group, display, lookup, gene_key, mut_key in _MATRIX_COLS:
            presence[sid][display] = _mut_present(res, lookup, gene_key, mut_key)

    active_muts = sorted(_col_order.keys(), key=lambda m: _col_order[m])

    with st.expander("Resistance mutation matrix", expanded=True):

        samples = list(valid.keys())
        seg_seq = [(
            _col_group[m],
            _is_plasmid(_col_lookup[m]),
            m,
        ) for m in active_muts]

        grp_spans: list[tuple[str, int]] = [
            (g, sum(1 for _ in it))
            for g, it in _groupby(seg_seq, key=lambda x: x[0])
        ]

        type_spans: list[tuple[str, bool, int]] = [
            (g, p, sum(1 for _ in it))
            for (g, p), it in _groupby(seg_seq, key=lambda x: (x[0], x[1]))
        ]

        sample_th = (
            '<th rowspan="3" style="background:#f1f5f9;color:#475569;font-size:0.7rem;'
            'padding:6px 12px;border:1px solid #d1d5db;text-align:left;'
            'font-weight:600;vertical-align:bottom;white-space:nowrap;">Sample</th>'
        )
        grp_cells = ""
        for grp, span in grp_spans:
            col = _GROUP_COLORS.get(grp, "#6b7280")
            lbl = _GROUP_LABELS.get(grp, grp)
            grp_cells += (
                f'<th colspan="{span}" style="background:{col};color:#fff;'
                f'font-size:0.65rem;font-weight:700;text-align:center;'
                f'padding:4px 6px;border:1px solid #fff;white-space:nowrap;">'
                f'{lbl}</th>'
            )

        type_cells = ""
        for grp, is_plasm, span in type_spans:
            col = _GROUP_COLORS.get(grp, "#6b7280")
            lbl = "Plasmid" if is_plasm else "Chromosomal"
            type_cells += (
                f'<th colspan="{span}" style="background:{col}18;color:{col};'
                f'font-size:0.6rem;font-weight:700;text-align:center;'
                f'padding:3px 4px;border:1px solid #d1d5db;'
                f'border-top:2px solid {col};white-space:nowrap;">'
                f'{lbl}</th>'
            )

        mut_cells = ""
        for m in active_muts:
            col = _GROUP_COLORS.get(_col_group[m], "#6b7280")
            mut_cells += (
                f'<th style="background:#f8fafc;border:1px solid #d1d5db;'
                f'padding:4px 2px;vertical-align:bottom;height:90px;width:28px;">'
                f'<div style="writing-mode:vertical-rl;transform:rotate(180deg);'
                f'color:#334155;font-size:0.62rem;font-weight:600;white-space:nowrap;'
                f'padding-bottom:4px;">{m}</div>'
                f'</th>'
            )

        data_rows = ""
        for sid in samples:
            cells = (
                f'<td style="color:#1e293b;font-size:0.75rem;font-weight:600;'
                f'padding:5px 12px;border:1px solid #d1d5db;white-space:nowrap;'
                f'background:#f8fafc;">{sid}</td>'
            )
            for m in active_muts:
                col = _GROUP_COLORS.get(_col_group[m], "#6b7280")
                if presence[sid].get(m):
                    cells += (
                        f'<td style="text-align:center;border:1px solid #e2e8f0;'
                        f'padding:5px;background:#fff;" title="{m}">'
                        f'<span style="color:{col};font-size:1rem;line-height:1;">&#9679;</span>'
                        f'</td>'
                    )
                else:
                    cells += '<td style="background:#fff;border:1px solid #e2e8f0;padding:5px;"></td>'
            data_rows += f"<tr>{cells}</tr>"

        table_html = f"""
<div style="overflow-x:auto;margin-bottom:0.5rem;">
<table style="border-collapse:collapse;font-family:sans-serif;background:#fff;min-width:100%;">
  <thead>
    <tr>{sample_th}{grp_cells}</tr>
    <tr>{type_cells}</tr>
    <tr>{mut_cells}</tr>
  </thead>
  <tbody>{data_rows}</tbody>
</table>
</div>
<div style="font-size:0.68rem;color:#64748b;margin-top:4px;">
  &#9679; = mutation/gene present &nbsp;·&nbsp; empty cell = absent
</div>"""
        st.markdown(table_html, unsafe_allow_html=True)

        csv_rows = []
        for sid in samples:
            row: dict = {"Sample": sid}
            for m in active_muts:
                row[m] = "●" if presence[sid].get(m) else ""
            csv_rows.append(row)
        st.download_button(
            "⬇ Download matrix (CSV)",
            data=pd.DataFrame(csv_rows).to_csv(index=False).encode(),
            file_name="mutation_matrix.csv",
            mime="text/csv",
            key=f"mut_matrix_csv_{id(results_dict)}",
        )


_ESSENTIAL_GENE_INFO = {
    "pilT":  ("Motility (type IV pilus)", "Twitching motility — retraction ATPase"),
    "pilT2": ("Motility (type IV pilus)", "Twitching motility — retraction ATPase paralogue"),
    "ftsZ":  ("Cell division", "Septum formation — tubulin homologue"),
    "recA":  ("DNA repair / recombination", "Homologous recombination, SOS response, natural transformation"),
    "comA":  ("Natural transformation (competence)", "DNA uptake across the outer membrane"),
    "tonB":  ("Iron acquisition", "Energy transducer — powers TonB-dependent outer membrane iron transporters"),
    "fur":   ("Iron acquisition / stress response", "Master regulator of iron-uptake and oxidative-stress genes"),
    "rpoH":  ("Stress response", "Sigma-32 factor — activates heat-shock gene expression"),
    "tbpB":  ("Iron acquisition", "Transferrin-binding protein — extracts iron from host transferrin"),
}


def _render_essential_gene_markers(results_dict: dict) -> None:
    valid = {sid: r for sid, r in results_dict.items() if isinstance(r, AMRResult)}
    if not valid:
        return

    rows = []
    for sid, res in valid.items():
        nonsyn = res.essential_gene_mutations or {}
        syn    = res.essential_gene_synonymous or {}
        for gene in sorted(_ESSENTIAL_GENE_INFO):
            process, product = _ESSENTIAL_GENE_INFO[gene]
            muts     = nonsyn.get(gene) or []
            syn_muts = syn.get(gene) or []
            if not muts and not syn_muts:
                continue
            rows.append({
                "Sample":      sid,
                "Gene":        gene,
                "Process":     process,
                "Product":     product,
                "Mutations":   ", ".join(muts) if muts else "—",
                "Synonymous":  len(syn_muts),
            })

    with st.expander("Other Gene Variants (non-AMR)", expanded=False):
        st.caption(
            "Mutations in genes with no known role in antimicrobial resistance — motility, cell "
            "division, DNA repair, natural transformation, iron acquisition, and stress response. "
            " It is shown for genomic context only: a variant here can influence "
            "fitness, transformation efficiency, or other processes unrelated to drug resistance."
        )
        if not rows:
            st.info("No variants detected in the tracked essential genes for these samples.")
            return
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
        st.download_button(
            "⬇ Download essential gene variants (CSV)",
            data=pd.DataFrame(rows).to_csv(index=False).encode(),
            file_name="essential_gene_variants.csv",
            mime="text/csv",
            key=f"essential_genes_csv_{id(results_dict)}",
        )


def _render_cgmlst_clustering(cgmlst_result: dict | None, key_prefix: str = "cg") -> None:
    if not cgmlst_result:
        return
    st.divider()
    st.subheader("cgMLST Clustering")
    st.caption(
        "Core genome clustering, source: [pubmlst.org](https://pubmlst.org/): "
        "*N. gonorrhoeae* cgMLST v1.0"
    )
    st.caption(
        f"Core genome: {cgmlst_result.get('core_loci_count', '?')} loci (≥95% presence). "
        "Single-linkage clustering at a fixed threshold, "
        "enabling longitudinal comparison."
    )

    _cg_stats = cgmlst_result.get("sample_stats", {})

    if not _cg_stats:
        return

    st.caption(
        "**400 allele difference threshold (Ng_cgc_400)**, stable genogroup "
    )
    _gg_rows = []
    for _sid, _st in _cg_stats.items():
        _gg_rows.append({
            "Sample":            _sid,
            "Genogroup/Cluster": _st.get("genogroup"),
            "% assigned":        f"{_st.get('pct_assigned', 0):.1f}%",
        })
    _gg_rows.sort(key=lambda r: (
        r["Genogroup/Cluster"] if r["Genogroup/Cluster"] is not None else float("inf"),
        r["Sample"],
    ))
    for _row in _gg_rows:
        _row["Genogroup/Cluster"] = _row["Genogroup/Cluster"] if _row["Genogroup/Cluster"] is not None else "—"
    st.dataframe(pd.DataFrame(_gg_rows), use_container_width=True, hide_index=True)

    _cg_pairs = cgmlst_result.get("pairwise", [])
    if _cg_pairs:
        with st.expander("Pairwise allele differences"):
            _pair_rows = []
            for _p in sorted(_cg_pairs, key=lambda x: x["allele_diff"]):
                _pair_rows.append({
                    "Sample A":          _p["sample_a"],
                    "Sample B":          _p["sample_b"],
                    "Allele differences": _p["allele_diff"],
                    "Same genogroup (≤400)": "Yes" if _p["same_genogroup"] else "No",
                })
            st.dataframe(pd.DataFrame(_pair_rows), use_container_width=True, hide_index=True)

    _cg_matrix_csv = cgmlst_result.get("matrix_csv")
    if _cg_matrix_csv and Path(_cg_matrix_csv).exists():
        st.download_button(
            "⬇ cgMLST distance matrix (CSV)",
            data=Path(_cg_matrix_csv).read_bytes(),
            file_name="cgmlst_distances.csv",
            mime="text/csv",
            key=f"{key_prefix}_matrix_dl",
        )


def _render_cohort_amr_profile(amr_samples: dict, project_name: str = "") -> None:
    import plotly.graph_objects as go
    import datetime as _dt

    _COHORT_MIN = 5
    valid = {sid: r for sid, r in amr_samples.items() if isinstance(r, AMRResult)}
    n = len(valid)

    if n < _COHORT_MIN:
        st.info(
            f"Cohort visualisation will appear when ≥{_COHORT_MIN} samples have been analysed "
            f"in this project ({n}/{_COHORT_MIN} ready)."
        )
        return

    active_cols = [
        c for c in _MATRIX_COLS
        if any(_mut_present(r, c[2], c[3], c[4]) for r in valid.values())
    ]
    if not active_cols:
        st.caption("No resistance determinants detected across the cohort.")
        return

    col_labels = [c[1] for c in active_cols]
    col_groups = [c[0] for c in active_cols]

    presence = {
        sid: {c[1]: _mut_present(r, c[2], c[3], c[4]) for c in active_cols}
        for sid, r in valid.items()
    }
    totals = {sid: sum(1 for v in presence[sid].values() if v) for sid in valid}
    sids = sorted(valid, key=lambda s: totals[s], reverse=True)

    _col_pheno: dict[str, str] = {}
    for c in active_cols:
        _, display, lookup, gene_key, mut_key = c
        if lookup in ("chrom", "promoter") and gene_key and mut_key:
            _col_pheno[display] = _CDC_RULES.get(f"{gene_key}_{mut_key}", c[0])
        elif lookup == "plasmid" and gene_key:
            _col_pheno[display] = _PLASMID_CDC_RULES.get(gene_key, c[0])
        elif lookup == "mosaic":
            _col_pheno[display] = "penA mosaic suspected (sequence identity, independent of the mutation-count check)"
        elif lookup == "truncation":
            _col_pheno[display] = "efflux_pump_overexpression"
        else:
            _col_pheno[display] = c[0]

    z_vals: list[list[int]] = []
    hover_txt: list[list[str]] = []
    for sid in sids:
        z_row, h_row = [], []
        for cl in col_labels:
            det = presence[sid].get(cl, False)
            pheno = _col_pheno.get(cl, "")
            z_row.append(1 if det else 0)
            h_row.append(
                f"<b>{sid}</b><br>{cl}<br><i>{pheno}</i>"
                if det else f"<b>{sid}</b><br>{cl}<br><i>absent</i>"
            )
        z_vals.append(z_row)
        hover_txt.append(h_row)

    st.divider()
    st.markdown(f"""
<div style="margin:1.2rem 0 0.8rem 0;padding-bottom:0.55rem;border-bottom:2px solid #e2e8f0;">
  <div style="font-size:1.05rem;font-weight:700;color:#0f172a;">Cohort AMR Profile: {n} samples</div>
</div>
""", unsafe_allow_html=True)

    n_s = len(sids)
    _row_h_px = max(22, min(48, 480 // max(1, n_s)))
    _fig_h = max(220, n_s * _row_h_px + 160)

    fig = go.Figure(go.Heatmap(
        z=z_vals,
        x=col_labels,
        y=sids,
        customdata=hover_txt,
        colorscale=[[0, "#1a2535"], [1, "#00c9b1"]],
        showscale=False,
        hovertemplate="%{customdata}<extra></extra>",
        zmin=0, zmax=1,
        xgap=2, ygap=2,
    ))

    for i, (g1, g2) in enumerate(zip(col_groups[:-1], col_groups[1:])):
        if g1 != g2:
            fig.add_vline(
                x=i + 0.5, line_width=2,
                line_color=_GROUP_COLORS.get(g2, "#4a5568"), opacity=0.6,
            )

    prev_grp, start_i = col_groups[0], 0
    for i, g in enumerate(col_groups[1:] + ["__end__"]):
        ri = i + 1
        if g != prev_grp:
            fig.add_annotation(
                x=(start_i + ri - 1) / 2.0, y=1.0, xref="x", yref="paper",
                yanchor="bottom", text=_GROUP_LABELS.get(prev_grp, prev_grp),
                showarrow=False, bgcolor=_GROUP_COLORS.get(prev_grp, "#4a5568"),
                borderpad=4, font=dict(color="white", size=8, family="monospace"), align="center",
            )
            prev_grp = g
            start_i = ri

    fig.update_layout(
        height=_fig_h,
        margin=dict(l=160, r=10, t=90, b=10),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        xaxis=dict(
            side="top", tickangle=-45,
            tickfont=dict(size=9, family="monospace", color="#8892a0"),
            showgrid=False,
        ),
        yaxis=dict(
            autorange="reversed",
            tickfont=dict(size=10, family="monospace", color="#8892a0"),
            showgrid=False,
        ),
        font=dict(color="#c8cdd5", family="monospace"),
    )

    st.plotly_chart(fig, use_container_width=True, key=f"cohort_profile_{id(amr_samples)}")
    st.caption(
        f"Source: European 2020 (IUSTI) AMR rules · {_dt.date.today().isoformat()} · {project_name or 'current project'}"
    )

    _csv_rows = []
    for sid in sids:
        row: dict = {"sample_id": sid}
        for cl in col_labels:
            row[cl] = 1 if presence[sid].get(cl, False) else 0
        row["total_determinants"] = totals[sid]
        pos_cats = sorted({
            _GROUP_LABELS.get(col_groups[i], col_groups[i])
            for i, cl in enumerate(col_labels)
            if presence[sid].get(cl, False)
        })
        row["resistance_categories"] = "; ".join(pos_cats) or "None"
        _csv_rows.append(row)

    st.download_button(
        "⬇ Download cohort AMR data (CSV)",
        data=pd.DataFrame(_csv_rows).to_csv(index=False).encode(),
        file_name=f"cohort_amr_{project_name or 'project'}.csv",
        mime="text/csv",
        key=f"cohort_csv_{id(amr_samples)}",
    )

    with st.expander("Cohort AMR summary charts"):
        cat_counts: dict[str, int] = {}
        for grp_code, grp_label in _GROUP_LABELS.items():
            idxs = [i for i, g in enumerate(col_groups) if g == grp_code]
            count = sum(
                1 for sid in sids
                if any(presence[sid].get(col_labels[i], False) for i in idxs)
            )
            if count > 0:
                cat_counts[grp_label] = count

        if cat_counts:
            _ca_df = pd.DataFrame([
                {"Category": k, "Count": v}
                for k, v in sorted(cat_counts.items(), key=lambda x: -x[1])
            ])
            fig_a = go.Figure(go.Bar(
                y=_ca_df["Category"],
                x=_ca_df["Count"],
                orientation="h",
                marker_color="#00c9b1",
                text=[f"{v}/{n}" for v in _ca_df["Count"]],
                textposition="outside",
                hovertemplate="%{y}: %{x} samples<extra></extra>",
            ))
            fig_a.update_layout(
                title=dict(text="Resistance category prevalence", font=dict(size=11, color="#c8cdd5", family="monospace")),
                height=max(220, len(cat_counts) * 36 + 100),
                margin=dict(l=160, r=80, t=50, b=20),
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)",
                xaxis=dict(
                    title=dict(text="Samples with ≥1 determinant", font=dict(color="#8892a0", size=9)),
                    range=[0, n * 1.2],
                    tickfont=dict(color="#8892a0", size=9),
                    showgrid=True, gridcolor="#1f2d3d",
                ),
                yaxis=dict(tickfont=dict(color="#8892a0", size=10, family="monospace"), showgrid=False),
                font=dict(color="#c8cdd5", family="monospace"),
            )
            st.plotly_chart(fig_a, use_container_width=True, key=f"cohort_cat_{id(amr_samples)}")

        def _bar_color(cnt: int) -> str:
            return "#dc2626" if cnt > 5 else ("#d97706" if cnt >= 3 else "#16a34a")

        fig_b = go.Figure(go.Bar(
            x=sids,
            y=[totals[sid] for sid in sids],
            marker_color=[_bar_color(totals[sid]) for sid in sids],
            customdata=[plural(totals[sid], "determinant") for sid in sids],
            hovertemplate="%{x}: %{y} %{customdata}<extra></extra>",
        ))
        fig_b.add_hline(
            y=3, line_dash="dash", line_color="#f97316", line_width=1.5,
            annotation_text="MDR threshold (≥3 classes)",
            annotation_position="top right",
            annotation_font=dict(color="#f97316", size=9),
        )
        fig_b.update_layout(
            title=dict(text="Per-sample resistance burden", font=dict(size=11, color="#c8cdd5", family="monospace")),
            height=max(220, 160 + n_s * 20),
            margin=dict(l=20, r=20, t=50, b=120),
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
            xaxis=dict(
                tickangle=-45,
                tickfont=dict(color="#8892a0", size=9, family="monospace"),
                showgrid=False,
            ),
            yaxis=dict(
                title=dict(text="# resistance determinants", font=dict(color="#8892a0", size=9)),
                tickfont=dict(color="#8892a0", size=9),
                showgrid=True, gridcolor="#1f2d3d",
            ),
            font=dict(color="#c8cdd5", family="monospace"),
        )
        st.plotly_chart(fig_b, use_container_width=True, key=f"cohort_burden_{id(amr_samples)}")


_RESISTANCE_CATEGORY_ORDER  = ["susceptible", "low_resistance", "moderate_resistance", "high_resistance", "MDR", "XDR"]
_RESISTANCE_CATEGORY_LABELS = {
    "susceptible":         "Susceptible",
    "low_resistance":      "Low resistance",
    "moderate_resistance": "Moderate resistance",
    "high_resistance":     "High resistance",
    "MDR":                 "MDR",
    "XDR":                 "XDR",
}
_RESISTANCE_CATEGORY_COLORS = {
    "Susceptible":         "#16a34a",
    "Low resistance":      "#84cc16",
    "Moderate resistance": "#facc15",
    "High resistance":     "#f97316",
    "MDR":                 "#dc2626",
    "XDR":                 "#7f1d1d",
    "Unknown":             "#cbd5e1",
}
_RESISTANCE_CATEGORY_LABEL_ORDER = [_RESISTANCE_CATEGORY_LABELS[c] for c in _RESISTANCE_CATEGORY_ORDER] + ["Unknown"]


def _sample_resistance_categories(project_name: str | None) -> dict[str, str]:
    """{sample_id: category label}."""
    if not project_name or db is None:
        return {}

    out: dict[str, str] = {}
    for sid, d in _load_project_cached(project_name).get("samples", {}).items():
        amr = d.get("amr")
        if amr is None:
            continue
        cat = amr.resistance_category or "susceptible"
        out[sid] = _RESISTANCE_CATEGORY_LABELS.get(cat, cat)
    return out


def _resistant_sample_ids(amr_samples: dict, group_code: str) -> set[str]:
    idxs = [c for c in _MATRIX_COLS if c[0] == group_code]
    return {
        sid for sid, r in amr_samples.items()
        if isinstance(r, AMRResult) and any(_mut_present(r, c[2], c[3], c[4]) for c in idxs)
    }


def _render_amr_tree_highlight(
    project_name: str | None, upload_names: list[str], key_prefix: str,
) -> dict[str, str] | None:
    """Antibiotic selectbox for tree highlighting. Returns {sample_id: "resistant"/
    "susceptible"/"no_data"} for uploaded samples under the chosen antibiotic, or None
    if "None" is selected."""
    if not upload_names:
        return None

    amr_samples: dict = {}
    if project_name and db is not None:
        amr_samples = {
            sid: d["amr"]
            for sid, d in _load_project_cached(project_name).get("samples", {}).items()
            if "amr" in d
        }

    groups_present = sorted({c[0] for c in _MATRIX_COLS}, key=lambda g: _GROUP_LABELS.get(g, g))
    options = ["None"] + [_GROUP_LABELS.get(g, g) for g in groups_present]
    choice = st.selectbox("Highlight resistance to", options, key=f"{key_prefix}_amr_highlight")
    if choice == "None":
        return None

    group_code = next(g for g in groups_present if _GROUP_LABELS.get(g, g) == choice)
    resistant = _resistant_sample_ids(amr_samples, group_code)

    result: dict[str, str] = {}
    for sid in upload_names:
        if sid not in amr_samples:
            result[sid] = "no_data"
        elif sid in resistant:
            result[sid] = "resistant"
        else:
            result[sid] = "susceptible"

    if not amr_samples:
        st.caption("No AMR data available yet for the uploaded samples.")
    elif "no_data" in result.values():
        st.caption("Samples without AMR data yet are shown greyed out.")
    return result


def render_metadata_charts(meta: dict, sample_ids: list[str], project_name: str | None = None) -> None:
    rows = [
        {"sample_id": sid, **meta[sid]}
        for sid in sample_ids
        if meta.get(sid) and any(meta[sid].values())
    ]
    if not rows:
        return

    cat_by_sid = _sample_resistance_categories(project_name)
    if not any(sid in cat_by_sid for sid in sample_ids):
        st.caption("Resistance breakdowns will appear here once AMR profiling has been run for these samples.")
        return

    df = pd.DataFrame(rows)
    df["Resistance category"] = df["sample_id"].map(cat_by_sid).fillna("Unknown")
    _cat_scale = alt.Scale(
        domain=_RESISTANCE_CATEGORY_LABEL_ORDER,
        range=[_RESISTANCE_CATEGORY_COLORS[c] for c in _RESISTANCE_CATEGORY_LABEL_ORDER],
    )

    date_col = df.get("collection_date")
    if date_col is not None and date_col.astype(bool).any():
        dt = df[df["collection_date"].astype(bool)].copy()
        dt["collection_date"] = pd.to_datetime(dt["collection_date"], errors="coerce", dayfirst=True)
        dt = dt.dropna(subset=["collection_date"])
        if not dt.empty:
            dt_counts = (
                dt.groupby([dt["collection_date"].dt.date, "Resistance category"])
                .size().reset_index(name="Samples")
            )
            dt_counts.columns = ["Date", "Resistance category", "Samples"]
            st.altair_chart(
                alt.Chart(dt_counts).mark_bar().encode(
                    x=alt.X("Date:T", title="Collection date"),
                    y=alt.Y("Samples:Q", title="Samples", axis=alt.Axis(tickMinStep=1)),
                    color=alt.Color("Resistance category:N", scale=_cat_scale,
                                     sort=_RESISTANCE_CATEGORY_LABEL_ORDER),
                    order=alt.Order("Resistance category:N", sort="ascending"),
                    tooltip=["Date", "Resistance category", "Samples"],
                ).properties(height=220, title="Samples by collection date, by resistance category"),
                use_container_width=True,
            )

    def _pct_resistant_chart(col: str, title: str):
        s = df.get(col)
        if s is None or not s.astype(bool).any():
            return None
        sub = df[s.astype(bool) & (df["Resistance category"] != "Unknown")]
        if sub.empty:
            return None
        grp = sub.groupby(col).agg(
            n=("Resistance category", "size"),
            pct_resistant=("Resistance category", lambda x: round((x != "Susceptible").mean() * 100, 1)),
        ).reset_index()
        grp.columns = [title, "n", "% resistant"]
        return alt.Chart(grp).mark_bar(color="#dc2626").encode(
            x=alt.X("% resistant:Q", title="% resistant", scale=alt.Scale(domain=[0, 100])),
            y=alt.Y(f"{title}:N", title=None, sort="-x"),
            tooltip=[title, "% resistant", alt.Tooltip("n:Q", title="Samples")],
        ).properties(height=160, title=f"% resistant by {title.lower()}")

    age_pct_chart = None
    age_s = df.get("age")
    if age_s is not None:
        ages_df = df[(df["Resistance category"] != "Unknown")].copy()
        ages_df["age_num"] = pd.to_numeric(ages_df.get("age"), errors="coerce")
        ages_df = ages_df.dropna(subset=["age_num"])
        if not ages_df.empty:
            ages_df["Age group"] = pd.cut(
                ages_df["age_num"], bins=[0, 20, 30, 40, 50, 200],
                labels=["<20", "20-29", "30-39", "40-49", "50+"], right=False,
            )
            age_grp = ages_df.groupby("Age group", observed=True).agg(
                n=("Resistance category", "size"),
                pct_resistant=("Resistance category", lambda x: round((x != "Susceptible").mean() * 100, 1)),
            ).reset_index()
            if not age_grp.empty:
                age_pct_chart = alt.Chart(age_grp).mark_bar(color="#dc2626").encode(
                    x=alt.X("Age group:N", title="Age", sort=["<20", "20-29", "30-39", "40-49", "50+"]),
                    y=alt.Y("pct_resistant:Q", title="% resistant", scale=alt.Scale(domain=[0, 100])),
                    tooltip=["Age group", alt.Tooltip("pct_resistant:Q", title="% resistant"),
                             alt.Tooltip("n:Q", title="Samples")],
                ).properties(height=160, title="% resistant by age")

    dist_charts = [
        c for c in (
            _pct_resistant_chart("anatomical_site", "Site"),
            _pct_resistant_chart("sex", "Sex"),
            age_pct_chart,
            _pct_resistant_chart("country", "Country"),
        ) if c is not None
    ]
    if dist_charts:
        for col, chart in zip(st.columns(len(dist_charts)), dist_charts):
            with col:
                st.altair_chart(chart, use_container_width=True)

def _inline_metadata_widget(
    sample_ids: list[str],
    project_name: str | None,
    key_prefix: str,
) -> None:
    """
    Optional metadata form shown after upload files in any module.
    Renders an expander.
    """
    if not sample_ids or not project_name or db is None:
        return

    existing = db.load_metadata(project_name) if db else {}

    def _parse_date(v):
        if not v:
            return None
        ts = pd.to_datetime(str(v), errors="coerce", dayfirst=True)
        return ts.date() if pd.notna(ts) else None

    with st.expander("Add metadata (optional)", expanded=False):

        blank = {f: "" for f in _META_COLS[1:]}
        rows = []

        for sid in sample_ids:
            row = {"sample_id": sid}
            row.update({
                **blank,
                **{k: (v or "") for k, v in existing.get(sid, {}).items()}
            })
            row["collection_date"] = _parse_date(
                existing.get(sid, {}).get("collection_date")
            )
            rows.append(row)

        edited = st.data_editor(
            pd.DataFrame(rows, columns=_META_COLS),
            column_config={
                "sample_id": st.column_config.TextColumn(
                    "Sample ID", disabled=True
                ),
                "collection_date": st.column_config.DateColumn(
                    "Collection date", format="YYYY-MM-DD"
                ),
                "anatomical_site": st.column_config.SelectboxColumn(
                    "Anatomical site", options=_SITE_OPTIONS
                ),
                "sex": st.column_config.SelectboxColumn(
                    "Sex", options=_SEX_OPTIONS
                ),
                "age": st.column_config.NumberColumn(
                    "Age", min_value=0, max_value=120, step=1
                ),
                "country": st.column_config.TextColumn("Country"),
                "city": st.column_config.TextColumn("City/Municipality"),
                "health_unit": st.column_config.TextColumn("Health unit"),
            },
            use_container_width=True,
            hide_index=True,
            num_rows="fixed",
            key=f"{key_prefix}_meta_editor",
        )

        if st.button("Save metadata", key=f"{key_prefix}_meta_save"):
            _save_df = edited.copy()

            def _serialise_date(v):
                if v is None or (isinstance(v, float) and pd.isna(v)):
                    return None

                if isinstance(v, str):
                    ts = pd.to_datetime(v, errors="coerce", dayfirst=False)
                    return ts.date().isoformat() if pd.notna(ts) else None

                return v.isoformat() if pd.notna(v) else None

            _save_df["collection_date"] = _save_df["collection_date"].apply(
                _serialise_date
            )

            db.init_project(project_name)
            db.save_metadata(
                project_name,
                _save_df.to_dict(orient="records")
            )

            _load_project_cached.clear()

            st.success(
                f"Metadata saved for {len(rows)} "
                f"{plural(len(rows), 'sample')}."
            )

            existing = db.load_metadata(project_name)

        render_metadata_charts(
            existing,
            sample_ids,
            project_name=project_name
        )

_PROJECT_RESULT_KEYS: tuple[str, ...] = (
    "full_pipeline_results",
    "amr_manual_results", "amr_manual_bytes", "amr_manual_errors",
    "assembly_results",
    "qc_manual_results", "fasta_qc_results",
    "last_tree_path", "last_cluster_report", "last_clusters",
    "last_matrix_path", "last_cgmlst_result", "last_upload_names",
    "phy_tree_path", "phy_run_id", "phy_cluster_report",
    "phy_matrix_path", "phy_upload_names", "phy_cgmlst_result",
    "phy_nj_tree_path", "phy_ml_tree_path",
    "_amr_loaded_project", "_asm_loaded_project", "_qc_loaded_project",
    "_amr_worker_pending_jobs", "_amr_worker_pending_project",
    "_asm_pending_jobs", "_asm_pending_project",
    "_qc_pending_jobs", "_qc_pending_project",
    "_fqc_pending_jobs", "_fqc_pending_project",
    "man_selected_sid", "pipe_selected_sid",
)

def _clear_project_results() -> None:
    for _k in _PROJECT_RESULT_KEYS:
        st.session_state.pop(_k, None)


def ensure_result_page_state() -> None:
    if "current_page" not in st.session_state:
        st.session_state["current_page"] = "Homepage"
    if "active_project" not in st.session_state:
        st.session_state["active_project"] = None
    if "full_pipeline_results" not in st.session_state:
        st.session_state["full_pipeline_results"] = None


PHENOTYPE_LABELS = {
    "ceftriaxone_reduced_susceptibility": ("Ceftriaxone — Reduced susceptibility",          "🔴"),
    "high_level_azithromycin_resistance": ("Azithromycin — High-level resistance",           "🔴"),
    "moderate_azithromycin_resistance":   ("Azithromycin — Moderate resistance",             "🟡"),
    "ciprofloxacin_resistant":            ("Ciprofloxacin — Resistant",                      "🔴"),
    "ciprofloxacin_intermediate_parC":    ("Ciprofloxacin — Intermediate (parC)",             "🟡"),
    "penicillin_resistance":              ("Penicillin — Beta-lactamase mediated (blaTEM-1)", "🔴"),
    "tetracycline_resistance":            ("Tetracycline — Plasmid-mediated (tet-M)",         "🟡"),
    "reduced_beta_lactam_susceptibility":  ("Beta-lactam — Reduced susceptibility (porB)",       "🟡"),
    "efflux_pump_overexpression":          ("Efflux pump overexpression (mtrR)",                 "🟡"),
    "tetracycline_chromosomal_resistance": ("Tetracycline — Chromosomal resistance (rpsJ V57M)", "🟡"),
    "reduced_penicillin_susceptibility":   ("Penicillin — Reduced susceptibility (ponA L421P)",  "🟡"),
    "zoliflodacin_reduced_susceptibility": ("Zoliflodacin — Reduced susceptibility (investigational)", "🟡"),
    "wildtype":                            ("No resistance detected — Wildtype",                 "🟢"),
}

_BADGE_PALETTE = {
    "🔴": ("#fef2f2", "#dc2626", "#fecaca"),
    "🟡": ("#fffbeb", "#d97706", "#fde68a"),
    "🟢": ("#f0fdf4", "#16a34a", "#bbf7d0"),
    "⚪": ("#f9fafb", "#6b7280", "#e5e7eb"),
}

_TREATMENT_DRUGS: list[tuple[str, str]] = [
    ("Ceftriaxone", "ceftriaxone"), ("Azithromycin", "azithromycin"),
    ("Gentamicin", "gentamicin"), ("Ciprofloxacin", "ciprofloxacin"),
    ("Penicillin", "penicillin"), ("Tetracycline", "tetracycline"),
    ("Doxycycline", "doxycycline"),
]


def _render_drug_cards(therapy: dict) -> None:
    avoid_low = [a.lower() for a in therapy.get("avoid", [])]
    rec_str   = " ".join(therapy.get("recommend", [])).lower()
    cards = ""
    for label, key in _TREATMENT_DRUGS:
        if key in avoid_low:
            bg, brd, txt, icon, sub = "#fef2f2", "#dc2626", "#991b1b", "✕", "Avoid"
        elif key in rec_str:
            bg, brd, txt, icon, sub = "#f0fdf4", "#16a34a", "#15803d", "✓", "Recommended"
        else:
            bg, brd, txt, icon, sub = "#f9fafb", "#d1d5db", "#6b7280", "—", "Not indicated"
        cards += (
            f'<div style="background:{bg};border:1.5px solid {brd};border-radius:10px;'
            f'padding:0.7rem 0.8rem;text-align:center;min-width:90px;flex:1;">'
            f'<div style="font-size:1.1rem;font-weight:700;color:{txt};">{icon}</div>'
            f'<div style="font-size:0.78rem;font-weight:700;color:{txt};margin:0.15rem 0;">{label}</div>'
            f'<div style="font-size:0.65rem;color:{txt};opacity:0.8;">{sub}</div>'
            f'</div>'
        )
    st.markdown(
        f'<div style="display:flex;gap:0.4rem;flex-wrap:wrap;margin-bottom:0.8rem;">{cards}</div>',
        unsafe_allow_html=True,
    )


def _render_drug_matrix(samples_therapy: dict[str, dict]) -> None:
    """One row per sample, one column per drug (recommended/avoid/not indicated), no expanding needed."""
    rows = []
    for sid, therapy in samples_therapy.items():
        therapy   = therapy or {}
        avoid_low = [a.lower() for a in therapy.get("avoid", [])]
        rec_str   = " ".join(therapy.get("recommend", [])).lower()
        row = {"Sample": sid}
        for label, key in _TREATMENT_DRUGS:
            if key in avoid_low:
                row[label] = "✕"
            elif key in rec_str:
                row[label] = "✓"
            else:
                row[label] = "—"
        rows.append(row)
    if not rows:
        return

    df = pd.DataFrame(rows)
    drug_cols = [label for label, _ in _TREATMENT_DRUGS]

    def _cell_style(v):
        if v == "✓":
            return "background-color:#f0fdf4;color:#16a34a;font-weight:700;"
        if v == "✕":
            return "background-color:#fef2f2;color:#dc2626;font-weight:700;"
        return "color:#9ca3af;"

    st.dataframe(
        df.style.applymap(_cell_style, subset=drug_cols),
        use_container_width=True,
        hide_index=True,
    )


def render_amr_interpretation(amr_result: AMRResult, sample_label: str = ""):
    cdc_list  = amr_result.cdc_phenotypes
    prob      = amr_result.failure_probability
    therapy   = amr_result.therapy
    chrom     = amr_result.chromosomal
    plasm     = amr_result.plasmid
    is_wt     = cdc_list == ["wildtype"]
    n_resist  = len([p for p in cdc_list if p != "wildtype"])

    #  Summary banner () gradient white to dark red)
    if is_wt:
        sbg, sc = "#f0fdf4", "#16a34a"
        sb = "#bbf7d0"
        slabel, sicon = "No resistance detected: Wildtype", "✅"
    else:
        sbg, sc = _resistance_color(prob)
        sb = sbg
        _mech = plural(n_resist, "mechanism")
        slabel = (f"High resistance risk: {n_resist} {_mech} detected"
                  if prob >= 0.4 else f"Resistance detected: {n_resist} {_mech}")
        sicon = "🔴" if prob >= 0.4 else "🟡"

    sample_html = (f"<div style='font-size:0.85rem;color:{sc};opacity:0.75;margin-bottom:0.3rem;'>"
                   f"Sample: <code>{sample_label}</code></div>") if sample_label else ""
    st.markdown(f"""
    <div style="border:1.5px solid {sb}; background:{sbg}; border-radius:12px;
                padding:1.1rem 1.4rem; margin-bottom:1rem;">
        {sample_html}
        <div style="font-size:1.25rem; font-weight:700; color:{sc};">{sicon} {slabel}</div>
        <div style="color:{sc}; opacity:0.7; font-size:0.82rem; margin-top:0.25rem;">
            European 2020 (IUSTI) &nbsp;·&nbsp; Genomic resistance score: {prob*100:.1f}%
        </div>
    </div>""", unsafe_allow_html=True)

    #  MDR/XDR strain alert 
    _who_matches = amr_result.who_matches
    if _who_matches:
        _n_strains = len(_who_matches)
        _strain_word = plural(_n_strains, "strain")
        _has_xdr = any(m["mdr_class"] == "XDR" for m in _who_matches)
        _alert_border = "#ef4444" if _has_xdr else "#f59e0b"
        _alert_bg     = "rgba(239,68,68,0.07)" if _has_xdr else "rgba(245,158,11,0.07)"
        _alert_label  = "XDR: Extensively Drug-Resistant" if _has_xdr else "MDR: Multidrug-Resistant"
        _alert_icon   = "🚨" if _has_xdr else "⚠️"
        _strain_tags  = "".join(
            f'<span style="display:inline-block;background:rgba(239,68,68,0.15);'
            f'border:1px solid rgba(239,68,68,0.4);border-radius:6px;'
            f'padding:0.1rem 0.5rem;font-size:0.75rem;font-weight:700;'
            f'color:#fca5a5;margin:0.15rem 0.2rem 0;white-space:nowrap;">'
            f'{m["strain"]} ({m["mdr_class"]})</span>'
            if m["mdr_class"] == "XDR" else
            f'<span style="display:inline-block;background:rgba(245,158,11,0.15);'
            f'border:1px solid rgba(245,158,11,0.35);border-radius:6px;'
            f'padding:0.1rem 0.5rem;font-size:0.75rem;font-weight:700;'
            f'color:#fcd34d;margin:0.15rem 0.2rem 0;white-space:nowrap;">'
            f'{m["strain"]} ({m["mdr_class"]})</span>'
            for m in _who_matches
        )
        _class_list = sorted({c for m in _who_matches for c in m["matched_classes"]})
        st.markdown(f"""
<div style="border:1.5px solid {_alert_border};background:{_alert_bg};border-radius:12px;
padding:1rem 1.4rem;margin-bottom:1rem;">
<div style="font-size:1rem;font-weight:700;color:{_alert_border};margin-bottom:0.4rem;">
{_alert_icon} Resistance profile consistent with WHO reference {_strain_word}: {_alert_label}
</div>
<div style="margin-bottom:0.5rem;">{_strain_tags}</div>
<div style="font-size:0.8rem;color:#94a3b8;line-height:1.6;">
Resistance classes matched: <strong style="color:#e2e8f0;">{", ".join(_class_list)}</strong><br>
This isolate meets the phenotypic criteria of the indicated WHO reference {_strain_word} from the
2012/2016 panels (Unemo et al. 2016, <em>J Antimicrob Chemother</em>).
Immediate clinical review is recommended.
</div>
</div>""", unsafe_allow_html=True)

    #  Colour scale legend 
    _scale_stops = [0.0, 0.15, 0.35, 0.55, 0.75, 1.0]
    _scale_html = "".join(
        f'<div style="flex:1;height:14px;background:{_resistance_color(p)[0]};'
        f'border-radius:{"4px 0 0 4px" if i==0 else ("0 4px 4px 0" if i==len(_scale_stops)-2 else "0")};"></div>'
        for i, p in enumerate(_scale_stops[:-1])
    )
    st.markdown(f"""
    <div style="margin-bottom:1rem;">
        <div style="font-size:0.75rem;color:#6b7280;margin-bottom:3px;">Resistance level</div>
        <div style="display:flex;height:14px;border-radius:4px;overflow:hidden;
                    border:1px solid #e5e7eb;">{_scale_html}</div>
        <div style="display:flex;justify-content:space-between;
                    font-size:0.68rem;color:#9ca3af;margin-top:2px;">
            <span>Susceptible</span><span>Low</span><span>Moderate</span>
            <span>High</span><span>Very high</span><span>Extreme</span>
        </div>
    </div>""", unsafe_allow_html=True)

    st.markdown("##### Resistance by agent")
    render_agent_resistance_table(amr_result)

    #  Mutations + plasmid
    col_chrom, col_plasm = st.columns([3, 2], gap="medium")

    with col_chrom:
        st.markdown("##### Chromosomal mutations")

        _t1_rows, _t2_rows, _t3_rows = [], [], []
        for g, muts in chrom.items():
            if g in _CHROM_GENE_SUPPRESS:
                continue
            gdisp = _CHROM_GENE_RENAME.get(g, g)
            for m in (muts or []):
                wkey = f"{g}_{m}"
                tier = _AMR_TIER.get(wkey, 2)
                row  = {"Gene": gdisp, "Mutation": m, "_gene": g, "_mut": m}
                if tier == 1:
                    _t1_rows.append(row)
                elif tier == 3:
                    _t3_rows.append(row)
                else:
                    _t2_rows.append(row)

        def _render_mut_table(rows: list) -> None:
            df = pd.DataFrame(rows)
            def _rs(row):
                impact = _mutation_impact(df.loc[row.name, "_gene"], df.loc[row.name, "_mut"])
                bg, txt = _impact_color(impact)
                return [f"background-color:{bg};color:{txt}"] * len(row)
            st.dataframe(df[["Gene", "Mutation"]].style.apply(_rs, axis=1),
                         use_container_width=True, hide_index=True)

        if not (_t1_rows or _t2_rows or _t3_rows):
            st.markdown("""<div style="background:#f0fdf4;border:1px solid #bbf7d0;
                border-radius:8px;padding:0.75rem 1rem;color:#15803d;font-size:0.9rem;">
                No chromosomal mutations detected</div>""", unsafe_allow_html=True)
        else:
            if _t1_rows:
                st.caption("Clinically actionable")
                _render_mut_table(_t1_rows)
            elif not _t2_rows:
                st.markdown("""<div style="background:#f0fdf4;border:1px solid #bbf7d0;
                    border-radius:8px;padding:0.75rem 1rem;color:#15803d;font-size:0.9rem;">
                    No Tier 1 mutations detected</div>""", unsafe_allow_html=True)

            if _t2_rows:
                with st.expander("Epidemiological markers", expanded=True):
                    _render_mut_table(_t2_rows)

            if _t3_rows:
                with st.expander("Advanced AMR details", expanded=False):
                    _render_mut_table(_t3_rows)

    with col_plasm:
        st.markdown("##### Plasmid determinants")
        for gene, status in plasm.items():
            present = status == "present"
            if present:
                bg, txt = _impact_color(_mutation_impact(gene))
                brd = bg
            else:
                bg, txt, brd = "#f9fafb", "#6b7280", "#e5e7eb"
            ico = "⊕ Present" if present else "○ Absent"
            st.markdown(f"""
            <div style="border:1px solid {brd};background:{bg};border-radius:8px;
                        padding:0.55rem 1rem;margin-bottom:0.45rem;
                        display:flex;justify-content:space-between;align-items:center;">
                <span style="font-weight:600;font-size:0.9rem;color:{txt};">{gene}</span>
                <span style="color:{txt};font-weight:600;font-size:0.88rem;">{ico}</span>
            </div>""", unsafe_allow_html=True)

    st.divider()

    #  Clinical recommendation
    avoid     = therapy.get("avoid", [])
    recommend = therapy.get("recommend", [])
    _avoid_low = [a.lower() for a in avoid]
    _rec_str   = " ".join(recommend).lower()

    # Prominent recommendation banner
    if recommend:
        _rec_display = "  ·  ".join(re.sub(r'\s+[\d(].*', '', r).strip().title() for r in recommend)
        _avoid_display = "  ·  ".join(re.sub(r'\s+[\d(].*', '', a).strip().title() for a in avoid) if avoid else None
        if prob >= 0.4:
            _rb_brd, _rb_accent = "#f59e0b", "#f59e0b"
            _rb_bg  = "linear-gradient(135deg,rgba(245,158,11,0.12),rgba(245,158,11,0.04))"
            _rb_lbl = "Treatment with caution, resistance detected"
        else:
            _rb_brd, _rb_accent = "#00c9b1", "#00c9b1"
            _rb_bg  = "linear-gradient(135deg,rgba(0,201,177,0.12),rgba(0,201,177,0.04))"
            _rb_lbl = "Recommended Treatment, European 2020 (IUSTI)"
        _avoid_html = (
            f'<div style="margin-top:0.6rem;font-size:0.82rem;color:#ef4444;">'
            f'<strong>Avoid:</strong> {_avoid_display}</div>'
        ) if _avoid_display else ""
        st.markdown(f"""
<div style="border:2px solid {_rb_brd};background:{_rb_bg};border-radius:14px;
        padding:1.2rem 1.5rem;margin:0.25rem 0 1rem 0;">
  <div style="font-size:0.72rem;font-weight:700;text-transform:uppercase;
          letter-spacing:0.1em;color:{_rb_accent};margin-bottom:0.5rem;">
💊 {_rb_lbl}
  </div>
  <div style="font-size:1.25rem;font-weight:800;color:#e2e8f0;line-height:1.35;">
{_rec_display}
  </div>
  <div style="font-size:0.78rem;color:#94a3b8;margin-top:0.35rem;">
Genomic resistance score: {prob*100:.1f}%
  </div>
  {_avoid_html}
</div>""", unsafe_allow_html=True)

        _render_drug_cards(therapy)

    _alternatives = therapy.get("alternatives", [])
    if _alternatives:
        _alt_rows = "".join(
            f'<div style="display:flex;gap:0.6rem;align-items:baseline;padding:0.2rem 0;">'
            f'<span style="font-size:0.72rem;font-weight:700;color:#64748b;min-width:185px;flex-shrink:0;">{a["condition"]}</span>'
            f'<span style="font-size:0.79rem;color:#94a3b8;">{a["regimen"]}</span>'
            f'</div>'
            for a in _alternatives
        )
        st.markdown(f"""
<div style="border:1px solid #334155;background:rgba(15,23,42,0.35);border-radius:10px;
     padding:0.75rem 1.2rem;margin:-0.2rem 0 0.8rem 0;">
  <div style="font-size:0.67rem;font-weight:700;text-transform:uppercase;
       letter-spacing:0.1em;color:#475569;margin-bottom:0.4rem;">Alternative Regimens</div>
  {_alt_rows}
</div>""", unsafe_allow_html=True)

    st.divider()

    #  CDC phenotype badges 
    st.markdown("##### Resistance Phenotypes")
    badges = ""
    for p in cdc_list:
        label, icon = PHENOTYPE_LABELS.get(p, (p, "⚪"))
        bg, txt, brd = _BADGE_PALETTE.get(icon, _BADGE_PALETTE["⚪"])
        badges += (f'<span style="display:inline-block;background:{bg};color:{txt};'
                   f'border:1px solid {brd};border-radius:20px;padding:0.28rem 0.85rem;'
                   f'margin:0.2rem;font-size:0.84rem;font-weight:500;">{icon} {label}</span>')
    st.markdown(f'<div style="margin-bottom:0.25rem;">{badges}</div>', unsafe_allow_html=True)

    st.divider()

    #  Genomic resistance score
    st.markdown("##### Genomic Resistance Score")
    bar_clr   = "#dc2626" if prob >= 0.4 else ("#f59e0b" if prob >= 0.15 else "#16a34a")
    risk_lbl  = "High" if prob >= 0.4 else ("Moderate" if prob >= 0.15 else "Low")

    g1, g2 = st.columns([1, 4])
    with g1:
        st.markdown(f"""
        <div style="text-align:center;padding:0.4rem 0;">
            <div style="font-size:2.1rem;font-weight:700;color:{bar_clr};">{prob*100:.1f}%</div>
            <div style="font-size:0.82rem;color:{bar_clr};font-weight:600;">{risk_lbl}</div>
        </div>""", unsafe_allow_html=True)
    with g2:
        st.markdown(f"""
        <div style="margin-top:1.1rem;">
            <div style="height:22px;background:#e5e7eb;border-radius:11px;overflow:hidden;">
                <div style="height:22px;width:{min(prob*100,100):.1f}%;
                            background:{bar_clr};border-radius:11px;"></div>
            </div>
            <div style="display:flex;justify-content:space-between;
                        font-size:0.72rem;color:#9ca3af;margin-top:0.2rem;">
                <span>0%</span><span>Low &lt;15%</span>
                <span>Moderate 15–40%</span><span>High ≥40%</span><span>100%</span>
            </div>
        </div>""", unsafe_allow_html=True)

    #  MLST 
    mlst_res = amr_result.mlst
    if mlst_res and "st" in mlst_res:
        st.divider()
        st.markdown("##### MLST — Neisseria Sequence Type")
        if mlst_res.get("error"):
            st.caption(f"MLST unavailable: {mlst_res['error']}")
        else:
            _st_val    = mlst_res.get("st") or "?"
            _novel     = mlst_res.get("novel", False)
            _inc       = mlst_res.get("incomplete", False)
            _alleles   = mlst_res.get("alleles", {})
            _st_known  = not _novel and not _inc and _st_val not in ("?", "new", "—", None)
            if _st_known:
                _ml_bg, _ml_brd, _st_color = "#eff6ff", "#bfdbfe", "#1d4ed8"
                _st_label = f"ST-{_st_val}"
            elif _novel:
                _ml_bg, _ml_brd, _st_color = "#fffbeb", "#fde68a", "#b45309"
                _st_label = "Novel"
            else:
                _ml_bg, _ml_brd, _st_color = "#fef2f2", "#fecaca", "#dc2626"
                _st_label = "?"
            _ml1, _ml2 = st.columns([1, 3])
            _ml1.markdown(f"""
            <div style="text-align:center;padding:0.6rem 0.5rem;
                        background:{_ml_bg};border:1.5px solid {_ml_brd};border-radius:10px;">
                <div style="font-size:1.6rem;font-weight:700;color:{_st_color};">{_st_label}</div>
                <div style="font-size:0.75rem;color:#6b7280;">Neisseria MLST</div>
                {f"<div style='font-size:0.72rem;color:#b45309;'>Novel {plural(sum(1 for v in _alleles.values() if v == 'new'), 'allele')}</div>" if _novel else ""}
                {"<div style='font-size:0.72rem;color:#dc2626;'>Incomplete profile</div>" if _inc and not _novel else ""}
            </div>""", unsafe_allow_html=True)
            if _alleles:
                _allele_df = pd.DataFrame(
                    [{"Locus": g, "Allele": a} for g, a in _alleles.items()]
                )
                _ml2.dataframe(_allele_df, use_container_width=True, hide_index=True)

    #  NG-STAR 
    ngstar_res = amr_result.ngstar
    if ngstar_res:
        st.divider()
        _NGSTAR_BASE = "https://pubmlst.org/organisms/neisseria-gonorrhoeae/ngstar"
        _NGSTAR_PROFILE = "https://pubmlst.org/bigsdb?db=pubmlst_neisseria_seqdef&page=profileInfo&scheme_id=67&profile_id={st}"
        st.markdown(
            "##### NG-STAR — Antimicrobial Resistance Sequence Type &nbsp;"
            f'<a href="{_NGSTAR_BASE}" target="_blank" rel="noopener" '
            'style="font-size:0.78rem;font-weight:400;color:#4b5563;'
            'border:1px solid #d1d5db;border-radius:4px;padding:2px 7px;'
            'text-decoration:none;vertical-align:middle;">'
            '🔗 PubMLST NG-STAR</a>',
            unsafe_allow_html=True,
        )
        if ngstar_res.get("error"):
            st.caption(f"NG-STAR unavailable: {ngstar_res['error']}")
        else:
            _ng_st      = ngstar_res.get("ST") or "?"
            _ng_alleles = ngstar_res.get("alleles", {})
            _ng_novel   = ngstar_res.get("novel", False)
            _ng_inc     = ngstar_res.get("incomplete", False)
            _st_known   = _ng_st not in ("?", "new", None) and not _ng_novel and not _ng_inc
            if _st_known:
                _ng_bg, _ng_brd, _ng_color = "#eff6ff", "#bfdbfe", "#1d4ed8"
                _ng_label = f"ST-{_ng_st}"
                _st_url   = _NGSTAR_PROFILE.format(st=_ng_st)
                _st_link  = f'<a href="{_st_url}" target="_blank" rel="noopener" style="font-size:0.72rem;color:#1d4ed8;text-decoration:underline;">Ver no PubMLST ↗</a>'
            elif _ng_novel:
                _ng_bg, _ng_brd, _ng_color = "#fffbeb", "#fde68a", "#b45309"
                _ng_label = "?"
                _st_link  = ""
            else:
                _ng_bg, _ng_brd, _ng_color = "#fef2f2", "#fecaca", "#dc2626"
                _ng_label = "?"
                _st_link  = ""
            _ngs1, _ngs2 = st.columns([1, 3])
            _ngs1.markdown(f"""
            <div style="text-align:center;padding:0.6rem 0.5rem;
                        background:{_ng_bg};border:1.5px solid {_ng_brd};border-radius:10px;">
                <div style="font-size:1.6rem;font-weight:700;color:{_ng_color};">{_ng_label}</div>
                <div style="font-size:0.75rem;color:#6b7280;margin-top:2px;">NG-STAR · PubMLST</div>
                {_st_link}
                {f"<div style='font-size:0.72rem;color:#b45309;margin-top:3px;'>Novel {plural(sum(1 for v in _ng_alleles.values() if v == 'new'), 'allele')}: ST undetermined</div>" if _ng_novel else ""}
                {"<div style='font-size:0.72rem;color:#dc2626;margin-top:3px;'>Incomplete profile</div>" if _ng_inc and not _ng_novel else ""}
            </div>""", unsafe_allow_html=True)
            if _ng_alleles:
                _gene_labels = {
                    "penA": "penA", "mtrR": "mtrR", "porB": "porB",
                    "ponA": "ponA", "gyrA": "gyrA", "parC": "parC", "23SrRNA": "23S rRNA",
                }
                _novel_genes = {g for g, a in _ng_alleles.items() if a in ("new", "error", "not_found", "no_ref")}
                _allele_rows = []
                for g, a in _ng_alleles.items():
                    _lbl = _gene_labels.get(g, g)
                    _flag = " ⚠ novel" if g in _novel_genes and a == "new" else (
                            " ✗ not found" if a in ("not_found", "no_ref") else "")
                    _allele_rows.append({"Gene": _lbl, "Allele": a, "Note": _flag.strip()})
                _df_ngs = pd.DataFrame(_allele_rows)
                if not _df_ngs["Note"].any():
                    _df_ngs = _df_ngs.drop(columns=["Note"])
                _ngs2.dataframe(_df_ngs, use_container_width=True, hide_index=True)

    #  NG-MAST
    ngmast_res = amr_result.ngmast
    if ngmast_res:
        st.divider()
        _NGMAST_BASE = "https://pubmlst.org/bigsdb?db=pubmlst_neisseria_seqdef&page=schemeInfo&scheme_id=71"
        _NGMAST_PROFILE = "https://pubmlst.org/bigsdb?db=pubmlst_neisseria_seqdef&page=profileInfo&scheme_id=71&profile_id={st}"
        st.markdown(
            "##### NG-MAST — Multi-Antigen Sequence Type &nbsp;"
        )
        if ngmast_res.get("error"):
            st.caption(f"NG-MAST unavailable: {ngmast_res['error']}")
        else:
            _nm_st      = ngmast_res.get("ST") or "?"
            _nm_alleles = ngmast_res.get("alleles", {})
            _nm_novel   = ngmast_res.get("novel", False)
            _nm_inc     = ngmast_res.get("incomplete", False)
            _nm_known   = _nm_st not in ("?", "new", None) and not _nm_novel and not _nm_inc
            if _nm_known:
                _nm_bg, _nm_brd, _nm_color = "#eff6ff", "#bfdbfe", "#1d4ed8"
                _nm_label = f"ST-{_nm_st}"
                _nm_url   = _NGMAST_PROFILE.format(st=_nm_st)
                _nm_link  = f'<a href="{_nm_url}" target="_blank" rel="noopener" style="font-size:0.72rem;color:#1d4ed8;text-decoration:underline;">Ver no PubMLST ↗</a>'
            elif _nm_novel:
                _nm_bg, _nm_brd, _nm_color = "#fffbeb", "#fde68a", "#b45309"
                _nm_label = "?"
                _nm_link  = ""
            else:
                _nm_bg, _nm_brd, _nm_color = "#fef2f2", "#fecaca", "#dc2626"
                _nm_label = "?"
                _nm_link  = ""
            _nm1, _nm2 = st.columns([1, 3])
            _nm1.markdown(f"""
            <div style="text-align:center;padding:0.6rem 0.5rem;
                        background:{_nm_bg};border:1.5px solid {_nm_brd};border-radius:10px;">
                <div style="font-size:1.6rem;font-weight:700;color:{_nm_color};">{_nm_label}</div>
                <div style="font-size:0.75rem;color:#6b7280;margin-top:2px;">NG-MAST · PubMLST</div>
                {_nm_link}
                {f"<div style='font-size:0.72rem;color:#b45309;margin-top:3px;'>Novel {plural(sum(1 for v in _nm_alleles.values() if v == 'new'), 'allele')}: ST undetermined</div>" if _nm_novel else ""}
                {"<div style='font-size:0.72rem;color:#dc2626;margin-top:3px;'>Incomplete profile</div>" if _nm_inc and not _nm_novel else ""}
            </div>""", unsafe_allow_html=True)
            if _nm_alleles:
                _nm_gene_labels = {"porB": "porB", "tbpB": "tbpB"}
                _nm_rows = []
                for g, a in _nm_alleles.items():
                    _lbl = _nm_gene_labels.get(g, g)
                    _flag = " ⚠ novel" if a == "new" else (
                            " ✗ not found" if a in ("not_found", "no_ref") else "")
                    _nm_rows.append({"Gene": _lbl, "Allele": a, "Note": _flag.strip()})
                _df_nm = pd.DataFrame(_nm_rows)
                if not _df_nm["Note"].any():
                    _df_nm = _df_nm.drop(columns=["Note"])
                _nm2.dataframe(_df_nm, use_container_width=True, hide_index=True)

    #  Mosaic penA Detection
    mosaic_res = amr_result.mosaic_pena
    if mosaic_res and "identity" in mosaic_res and mosaic_res["identity"] is not None:
        st.divider()
        st.markdown("##### Mosaic *penA* Detection")
        _mos_id    = mosaic_res.get("identity", "—")
        _mos_cov   = mosaic_res.get("coverage_pct", "—")
        _mos_cls   = mosaic_res.get("allele_class", "—")
        _mos_flag  = mosaic_res.get("mosaic_suspected", False)
        _mos_bg    = "#fef2f2" if _mos_flag else "#f0fdf4"
        _mos_brd   = "#dc2626" if _mos_flag else "#16a34a"
        _mos_txt   = "#991b1b" if _mos_flag else "#15803d"
        _mos_icon  = "⚠ Mosaic suspected" if _mos_flag else "✓ Non-mosaic"
        st.markdown(f"""
        <div style="border:1.5px solid {_mos_brd};background:{_mos_bg};border-radius:10px;
                    padding:0.9rem 1.2rem;margin-bottom:0.5rem;">
            <div style="font-weight:700;color:{_mos_txt};font-size:1rem;">{_mos_icon}</div>
            <div style="color:{_mos_txt};font-size:0.85rem;margin-top:0.3rem;">
                Allele class: <strong>{_mos_cls}</strong>
                &nbsp;·&nbsp; Identity to FA1090 penA: <strong>{_mos_id}%</strong>
                &nbsp;·&nbsp; Coverage: <strong>{_mos_cov}%</strong>
            </div>
            {"<div style='color:#991b1b;font-size:0.8rem;margin-top:0.4rem;'>Mosaic <em>penA</em> alleles (e.g. allele 60.001) confer high-level ceftriaxone resistance — MIC testing strongly recommended.</div>" if _mos_flag else ""}
        </div>""", unsafe_allow_html=True)
    elif mosaic_res.get("error"):
        pass  # silently skip if blastn not available

    syn_vars = amr_result.synonymous_variants

    _GENE_PRODUCT = {
        "penA":           "Penicillin-binding protein 2 (PBP2)",
        "ponA":           "Penicillin-binding protein 1A (PBP1A)",
        "porB":           "Outer membrane porin B",
        "mtrR":           "MtrCDE efflux pump repressor",
        "mtrC":           "MtrCDE efflux pump — MtrC subunit",
        "norM":           "NorM efflux pump (MATE family)",
        "macA":           "MacAB-TolC efflux pump — MacA subunit",
        "macB":           "MacAB-TolC efflux pump — MacB subunit",
        "gyrA":           "DNA gyrase subunit A (fluoroquinolone target)",
        "gyrB":           "DNA gyrase subunit B (fluoroquinolone target)",
        "parC":           "Topoisomerase IV subunit C (fluoroquinolone target)",
        "parE":           "Topoisomerase IV subunit E (fluoroquinolone target)",
        "23SrRNA":        "23S ribosomal RNA (azithromycin target)",
        "16SrRNA":        "16S ribosomal RNA (aminoglycoside target)",
        "rpsJ":           "Ribosomal protein S10 (tetracycline resistance)",
        "rpsE":           "Ribosomal protein S5 (spectinomycin resistance)",
        "rplD":           "Ribosomal protein L4 (macrolide resistance)",
        "rplV":           "Ribosomal protein L22 (macrolide resistance)",
        "folP":           "Dihydropteroate synthase (sulfonamide target)",
        "rpoB":           "RNA polymerase beta subunit (rifampicin target)",
        "rpoD":           "RNA polymerase sigma factor",
        "pilQ":           "Outer membrane secretin PilQ",
    }

    _AA_THREE = {
        "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
        "E": "Glu", "Q": "Gln", "G": "Gly", "H": "His", "I": "Ile",
        "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
        "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    }
    syn_rows = []
    for _sg, _smuts in syn_vars.items():
        for _entry in (_smuts or []):
            if isinstance(_entry, dict):
                _nuc  = _entry["nuc"]
                _aa   = _entry.get("aa", "?")
                _cpos = _entry.get("pos", 0)
                _m = re.match(r"([ACGT])(\d+)([ACGT])", _nuc, re.IGNORECASE)
                _aa3 = _AA_THREE.get(_aa, _aa)
                _notation = (
                    f"c.{_m.group(2)}{_m.group(1)}>{_m.group(3)} ({_aa3}{_cpos}{_aa3})"
                    if _m else _nuc
                )
            else:
                _nuc      = _entry
                _notation = _nuc
            syn_rows.append({
                "Gene":     _sg,
                "Product":  _GENE_PRODUCT.get(_sg, "—"),
                "Notation": _notation,
                "_nuc":     _nuc,
            })
    st.divider()
    st.markdown("""
<div style="background:linear-gradient(135deg,#1e293b 0%,#0f172a 100%);
        border-left:4px solid #6366f1;border-radius:10px;
        padding:0.9rem 1.2rem;margin-bottom:0.8rem;">
  <div style="font-size:0.65rem;font-weight:700;color:#a5b4fc;
          letter-spacing:0.12em;text-transform:uppercase;margin-bottom:0.25rem;">Genomic background</div>
  <div style="font-size:1.1rem;font-weight:800;color:#0f172a;">Synonymous variants</div>
  <div style="font-size:0.8rem;color:#94a3b8;margin-top:0.2rem;">Nucleotide changes that do not alter the amino acid sequence — useful for strain typing and phylogenetic background.</div>
</div>
""", unsafe_allow_html=True)
    if syn_rows:
        _panel_freq = load_synonymous_panel()
        _panel_meta = synonymous_panel_meta()
        _n_panel    = _panel_meta.get("n_genomes", 0)

        if _n_panel > 0:
            _show_all = st.toggle(
                "Show all synonymous variants (including rare/private)",
                value=False,
                key=f"syn_show_all_{sample_label}",
                help="By default, only variants seen in ≥5% of the reference panel are shown.",
            )
            _MIN_FREQ = 0.05
            _freq_filtered = (
                syn_rows if _show_all
                else [
                    r for r in syn_rows
                    if _panel_freq.get(r["Gene"], {}).get(r["_nuc"], 0.0) >= _MIN_FREQ
                ]
            )
        else:
            _freq_filtered = syn_rows

        _scol1, _scol2 = st.columns([2, 2])
        with _scol1:
            _all_genes = sorted({r["Gene"] for r in syn_rows})
            _gene_filter = st.multiselect(
                "Filter by gene",
                options=_all_genes,
                default=[],
                key=f"syn_gene_filter_{sample_label}",
                placeholder="All genes",
            )
        with _scol2:
            _search = st.text_input(
                "Search mutation",
                value="",
                key=f"syn_search_{sample_label}",
                placeholder="e.g. c.273C>T",
            )

        filtered = [
            r for r in _freq_filtered
            if (not _gene_filter or r["Gene"] in _gene_filter)
            and (not _search or _search.upper() in r["Notation"].upper()
                 or _search.upper() in r["_nuc"].upper())
        ]

        if _n_panel > 0:
            for r in filtered:
                _fv = _panel_freq.get(r["Gene"], {}).get(r["_nuc"])
                r["Panel freq (%)"] = f"{_fv * 100:.1f}" if _fv is not None else "—"
            col_order = ["Gene", "Product", "Notation", "Panel freq (%)"]
            st.caption(f"Panel frequency from {_n_panel} reference {plural(_n_panel, 'genome')}. Variants absent from the panel are shown as —.")
        else:
            col_order = ["Gene", "Product", "Notation"]
            st.caption(
                "No reference panel loaded. Run `python backend/build_synonymous_panel.py "
                "--genomes /path/to/reference/genomes/` to populate frequency data."
            )

        if filtered:
            st.dataframe(
                pd.DataFrame(filtered)[col_order],
                use_container_width=True,
                hide_index=True,
            )
            st.caption(f"{len(filtered)} synonymous {plural(len(filtered), 'variant')} shown · {len(syn_rows)} total detected")
        else:
            st.info("No variants match the current filter.")
    else:
        st.caption("No synonymous variants detected in the screened genes for this sample.")

#  Overview panel: metrics + score chart

__all__ = [n for n in dir() if not n.startswith("__")]
