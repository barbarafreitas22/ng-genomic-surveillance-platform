from datetime import datetime
from pathlib import Path


def _plural(n: int, word: str, plural_form: str | None = None) -> str:
    return word if n == 1 else (plural_form or f"{word}s")


_PHENOTYPE_LABEL = {
    "ceftriaxone_reduced_susceptibility":  "Ceftriaxone reduced susceptibility",
    "high_level_azithromycin_resistance":  "Azithromycin high-level resistance",
    "moderate_azithromycin_resistance":    "Azithromycin moderate resistance",
    "ciprofloxacin_resistant":             "Ciprofloxacin resistance",
    "ciprofloxacin_intermediate_parC":     "Ciprofloxacin intermediate",
    "penicillin_resistance":               "Penicillin resistance",
    "tetracycline_resistance":             "Tetracycline resistance",
    "reduced_beta_lactam_susceptibility":  "Beta-lactam reduced susceptibility",
    "tetracycline_chromosomal_resistance": "Tetracycline chromosomal resistance",
    "reduced_penicillin_susceptibility":   "Penicillin reduced susceptibility",
    "efflux_pump_overexpression":          "MtrCDE efflux: reduced susceptibility (penicillin, tetracycline, azithromycin)",
    "norM_efflux_upregulation":            "NorM efflux: reduced fluoroquinolone susceptibility",
    "efflux_tetracycline_contribution":    None,
    "penA_allele_likely_mosaic":           "Mosaic penA pattern suspected (≥3 associated mutations)",
    "fluoroquinolone_minor_parE":          None,
    "zoliflodacin_reduced_susceptibility": "Zoliflodacin reduced susceptibility (investigational)",
}

def _fmt_phenotypes(cdc: list, sep: str = "; ") -> str:
    labels = []
    for p in cdc:
        if p == "wildtype":
            continue
        if p in _PHENOTYPE_LABEL:
            lbl = _PHENOTYPE_LABEL[p]
        else:
            lbl = p.replace("_", " ").capitalize()
        if lbl is not None:
            labels.append(lbl)
    return sep.join(labels) or "Wildtype"


_RESISTANCE_COLORS = {
    "susceptible":         ("#f0fdf4", "#166534"),
    "low_resistance":      ("#fefce8", "#854d0e"),
    "moderate_resistance": ("#fff7ed", "#9a3412"),
    "high_resistance":     ("#fef2f2", "#991b1b"),
    "MDR":                 ("#fdf2f8", "#701a75"),
}

_QC_ICONS   = {"pass": "✅", "warn": "⚠️", "fail": "❌"}
_QC_COLORS  = {
    "pass": ("#f0fdf4", "#166534"),
    "warn": ("#fffbeb", "#92400e"),
    "fail": ("#fef2f2", "#991b1b"),
}

_ASM_THR_DEFS = [
    ("min_genome_fraction", "Core genes %",      "completeness", "%",  True),
    ("min_n50",             "N50",               "n50",          " bp", True),
    ("max_contigs",         "Contigs",           "n_contigs",    "",   False),
    ("min_total_len",       "Min total length",  "total_len",    " bp", True),
    ("max_total_len",       "Max total length",  "total_len",    " bp", False),
]


# ── HTML helpers ──────────────────────────────────────────────────────────────

def _badge(text: str, bg: str, fg: str) -> str:
    return (
        f'<span style="display:inline-block;padding:2px 8px;border-radius:999px;'
        f'font-size:11px;font-weight:600;background:{bg};color:{fg}">{text}</span>'
    )


def _resistance_badge(category: str) -> str:
    bg, fg = _RESISTANCE_COLORS.get(category, ("#f3f4f6", "#374151"))
    label = category.replace("_", " ").title() if category else "—"
    return _badge(label, bg, fg)


def _qc_badge(status: str) -> str:
    if not status:
        return "—"
    icon = _QC_ICONS.get(status, "")
    bg, fg = _QC_COLORS.get(status, ("#f3f4f6", "#374151"))
    return _badge(f"{icon} {status.upper()}", bg, fg)


def _pass_fail(passed: bool) -> str:
    return _badge("✅ Pass", "#f0fdf4", "#166534") if passed else _badge("❌ Fail", "#fef2f2", "#991b1b")


def _table(headers: list[str], rows: list[list], striped: bool = True) -> str:
    if not rows:
        return "<p style='font-size:12px;color:#6b7280;'>No data available.</p>"
    ths = "".join(f"<th>{h}</th>" for h in headers)
    trs = ""
    for i, row in enumerate(rows):
        bg = ' style="background:#f9fafb"' if striped and i % 2 == 1 else ""
        cells = "".join(f"<td>{cell}</td>" for cell in row)
        trs += f"<tr{bg}>{cells}</tr>"
    return (
        f'<table style="width:100%;border-collapse:collapse;font-size:12px;margin-bottom:16px">'
        f'<thead><tr style="background:#f3f4f6">{ths}</tr></thead>'
        f'<tbody>{trs}</tbody></table>'
    )


def _section(title: str, body: str, note: str = "") -> str:
    note_html = f'<p style="font-size:11px;color:#6b7280;margin-bottom:8px">{note}</p>' if note else ""
    return (
        f'<h2 style="font-size:15px;margin:28px 0 8px;color:#1e3a5f;'
        f'border-bottom:2px solid #e5e7eb;padding-bottom:4px">{title}</h2>'
        f'{note_html}{body}'
    )


def _css() -> str:
    return """
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
    font-size: 14px; color: #111827; background: #fff;
    max-width: 1100px; margin: 0 auto; padding: 24px 32px;
}
h1 { font-size: 20px; margin-bottom: 4px; color: #1e3a5f; }
th { text-align: left; padding: 6px 10px;
     font-size: 11px; text-transform: uppercase; letter-spacing: 0.05em; color: #374151; }
td { padding: 5px 10px; border-bottom: 1px solid #f3f4f6; vertical-align: top; }
tr:last-child td { border-bottom: none; }
@media print {
    body { padding: 10px; }
    h2 { page-break-before: auto; }
}
"""


# ── Section builders ──────────────────────────────────────────────────────────

def _build_overview(samples: dict) -> str:
    from backend.assembly import parse_contigs_stats
    from backend.phylogeny.cgmlst import core_genome_gene_count

    headers = [
        "Sample", "Coverage", "≥10× breadth", "Species", "Assembly QC",
        "Total length", "N50", "Core genes %", "ST (MLST)",
        "Resistance category", "CDC phenotype",
    ]
    rows = []
    for sid, data in sorted(samples.items()):
        qr  = data.get("qc", {}) or {}
        ar  = data.get("assembly", {}) or {}
        mr  = data.get("amr", {}) or {}
        met = (qr.get("metrics") or {}) if not qr.get("error") else {}
        k2  = (qr.get("kraken2") or {}) if not qr.get("error") else {}

        cov   = met.get("mean_coverage")
        b10x  = met.get("pct_breadth_10x")
        sp    = k2.get("species_status", "")
        sp_icon = {"confirmed": "✅", "likely": "⚠️", "contaminated": "❌"}.get(sp, "")

        ctg_path = ar.get("contigs_path") if not ar.get("error") else None
        cst = None
        if ctg_path and Path(ctg_path).exists():
            try:
                core_genome_gene_count(ctg_path)
                cst = parse_contigs_stats(Path(ctg_path))
            except Exception:
                pass

        from .models import AMRResult as _AMRResult
        _mr_valid = isinstance(mr, _AMRResult)
        asm_qc = _qc_badge(cst.qc_status if cst is not None else None)
        mlst   = mr.mlst if _mr_valid else {}
        st_str = f"ST-{mlst['st']}" if mlst.get("st") and not mlst.get("error") else "—"
        res    = mr.resistance_category if _mr_valid else None
        cdc    = mr.cdc_phenotypes if _mr_valid else []
        pheno  = _fmt_phenotypes(cdc)

        rows.append([
            f"<strong>{sid}</strong>",
            f"{cov:.1f}×"  if cov  is not None else "—",
            f"{b10x:.1f}%" if b10x is not None else "—",
            f"{sp_icon} {sp.replace('_',' ').capitalize()}" if sp else "—",
            asm_qc,
            f"{cst.total_len/1e6:.2f} Mb" if cst is not None else "—",
            f"{cst.n50/1e3:.1f} kb"        if cst is not None else "—",
            f"{cst.completeness:.1f}%"      if cst is not None else "—",
            st_str,
            _resistance_badge(res) if res else "—",
            pheno,
        ])
    return _table(headers, rows)


def _build_assembly_thresholds(samples: dict) -> str:
    from backend.assembly import parse_contigs_stats, ASSEMBLY_QC_THRESHOLDS
    from backend.phylogeny.cgmlst import core_genome_gene_count

    rows = []
    for sid, data in sorted(samples.items()):
        ar = data.get("assembly", {}) or {}
        ctg_path = ar.get("contigs_path") if not ar.get("error") else None
        if not ctg_path or not Path(ctg_path).exists():
            continue
        try:
            core_genome_gene_count(ctg_path)
            cst = parse_contigs_stats(Path(ctg_path))
        except Exception:
            continue
        if cst is None:
            continue
        _cst_d = cst.to_dict()
        _pass_tier = ASSEMBLY_QC_THRESHOLDS.get("pass", {})
        for thr_key, label, stat_key, unit, is_min in _ASM_THR_DEFS:
            thr = ASSEMBLY_QC_THRESHOLDS.get(thr_key, _pass_tier.get(thr_key))
            val = _cst_d.get(stat_key)
            if thr is None or val is None:
                continue
            passed = (val >= thr) if is_min else (val <= thr)
            sign   = "≥" if is_min else "≤"
            rows.append([
                sid, label,
                f"{val:,.0f}{unit}",
                f"{sign}{thr:,.0f}{unit}",
                _pass_fail(passed),
            ])
    return _table(["Sample", "Metric", "Value", "Threshold", "Status"], rows)


def _build_amr_mutations(samples: dict) -> str:
    from .models import AMRResult as _AMRResult
    rows = []
    for sid, data in sorted(samples.items()):
        mr = data.get("amr") or {}
        if not isinstance(mr, _AMRResult):
            continue
        for gene, muts in mr.chromosomal.items():
            for mut in (muts or []):
                rows.append([sid, "Chromosomal", gene, mut])
        for gene, status in mr.plasmid.items():
            if status == "present":
                rows.append([sid, "Plasmid", gene, "present"])
    return _table(["Sample", "Type", "Gene", "Mutation/Status"], rows)


def _build_clinical(samples: dict) -> str:
    from .models import AMRResult as _AMRResult
    html = ""
    for sid, data in sorted(samples.items()):
        mr = data.get("amr") or {}
        if not isinstance(mr, _AMRResult):
            continue
        res_cat = mr.resistance_category or "susceptible"
        bg, fg  = _RESISTANCE_COLORS.get(res_cat, ("#f3f4f6", "#374151"))
        prob    = mr.failure_probability or 0.0
        ther    = mr.therapy or {}
        rec     = ", ".join(ther.get("recommend") or []) or "—"
        avoid   = ", ".join(ther.get("avoid") or []) or "None"
        cdc     = mr.cdc_phenotypes or []
        pheno   = _fmt_phenotypes(cdc)
        mlst    = mr.mlst or {}
        st_str  = f"ST-{mlst['st']}" if mlst.get("st") and not mlst.get("error") else "—"
        n_res   = mr.n_resistance_classes or 0

        html += (
            f'<div style="border:1px solid #e5e7eb;border-radius:8px;padding:12px 16px;'
            f'margin-bottom:8px;border-left:4px solid {fg};background:#fafafa">'
            f'<div style="font-weight:700;margin-bottom:4px">{sid} '
            f'<span style="font-weight:400;font-size:11px;color:#6b7280">'
            f'MLST: {st_str} &nbsp;·&nbsp; Score: {prob*100:.1f}% &nbsp;·&nbsp; '
            f'{n_res} resistance class(es)</span></div>'
            f'<div style="font-size:12px;margin-bottom:4px"><strong>Phenotype:</strong> {pheno}</div>'
            f'<div style="display:flex;gap:24px;font-size:12px">'
            f'<span><strong style="color:#dc2626">⊗ Avoid:</strong> {avoid}</span>'
            f'<span><strong style="color:#16a34a">✓ Recommend:</strong> {rec}</span>'
            f'</div></div>'
        )
    return html or "<p style='font-size:12px;color:#6b7280'>No AMR data available.</p>"


def _build_phylogeny(phylogeny: dict | None) -> str:
    if not phylogeny or not phylogeny.get("cluster_report"):
        return "<p style='font-size:12px;color:#6b7280'>No phylogeny data available.</p>"

    cr = phylogeny["cluster_report"]
    has_pp = any("poppunk_cluster" in v for v in cr.values())
    headers = ["Sample", "Type", "Genogroup (≤2000 SNP)"]
    if has_pp:
        headers.append("PopPUNK Cluster")
    headers += ["Nearest neighbour", "Distance (SNPs)", "Distance (%)"]

    rows = []
    for sname, info in sorted(cr.items()):
        stype = "Uploaded" if info.get("is_uploaded") else "Reference"
        row = [
            f"<strong>{sname}</strong>" if info.get("is_uploaded") else sname,
            stype,
            str(info.get("cluster", "—")),
        ]
        if has_pp:
            row.append(info.get("poppunk_cluster", "—"))
        row += [
            info.get("nearest_neighbour", "—"),
            f"{info.get('nn_distance', 0):.0f}",
            f"{info.get('nn_distance_pct', 0):.4f}%",
        ]
        rows.append(row)
    return _table(headers, rows)


# ── Public API ────────────────────────────────────────────────────────────────

def generate_html_report(
    project_name: str,
    samples: dict,
    phylogeny: dict | None,
) -> str:
    """
    Generate a self-contained HTML report for a surveillance project.
    samples:   db.load_project()["samples"]
    phylogeny: db.load_project()["phylogeny"]  (may be None)
    """
    now      = datetime.now().strftime("%Y-%m-%d %H:%M")
    n_samp   = len(samples)
    from .models import AMRResult as _AMRResult
    n_amr    = sum(1 for d in samples.values() if isinstance(d.get("amr"), _AMRResult))
    n_asm    = sum(1 for d in samples.values() if d.get("assembly") and not d["assembly"].get("error"))

    overview_html   = _build_overview(samples)
    asm_thr_html    = _build_assembly_thresholds(samples)
    mut_html        = _build_amr_mutations(samples)
    clinical_html   = _build_clinical(samples)
    phylogeny_html  = _build_phylogeny(phylogeny)

    body = f"""
{_section("1. Genomic Overview", overview_html)}

{_section(
    "2. Assembly Quality — CDC AR Lab Network EQA Thresholds",
    asm_thr_html,
    note="Pass/caution tiers for <em>N. gonorrhoeae</em> (≥2 of 3 metrics): contigs ≤150/180, N50 &gt;30/20 kb, total length 2.0–2.2/1.8–2.2 Mb. Independent hard-fail checks: GC% 50–56%, core genome completeness ≥90% (BLASTN vs curated core gene set)."
)}

{_section("3. AMR Resistance Mutations", mut_html)}

{_section(
    "4. Clinical Interpretation",
    clinical_html,
    note="Based on European 2020 (IUSTI) treatment guidelines. Score = ordinal resistance severity index (0.0–1.0), not a calibrated probability."
)}

{_section(
    "5. Phylogeny — Cluster Report",
    phylogeny_html,
    note="Genogroup: cgMLST allele-based single-linkage clustering at ≤400 allele differences (Harrison et al. 2020), a fixed, dataset-independent threshold."
)}
"""

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>NG Surveillance Report: {project_name}</title>
<style>{_css()}</style>
</head>
<body>

<h1>&#x1F9EC; <em>Neisseria gonorrhoeae</em> Genomic Surveillance Report</h1>
<p style="font-size:12px;color:#6b7280;margin-bottom:4px">
  Project: <strong>{project_name}</strong> &nbsp;·&nbsp;
  {n_samp} {_plural(n_samp, 'sample')} &nbsp;·&nbsp;
  {n_asm} assembled &nbsp;·&nbsp;
  {n_amr} with AMR profiling &nbsp;·&nbsp;
  Generated: {now}
</p>
<p style="font-size:11px;color:#9ca3af;margin-bottom:24px">
  NG Genomic Surveillance Platform
</p>

{body}

<div style="margin-top:32px;padding-top:12px;border-top:1px solid #e5e7eb;
            font-size:11px;color:#9ca3af">
  Generated by the <em>N. gonorrhoeae</em> Genomic Surveillance Platform &nbsp;·&nbsp; {now}
</div>

</body>
</html>"""
