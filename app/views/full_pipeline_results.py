from __future__ import annotations
import json
import os
import re
import tempfile
import time
import base64
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

import altair as alt
import pandas as pd
import streamlit as st

from views.shared import *
from views.shared import (
    _load_project_cached, _resistance_color,
    _render_cohort_amr_profile, _render_genome_overview,
    _render_essential_gene_markers,
)


def _render_d3_tree(tree_path: str, clusters: dict, upload_names: list, legend_title: str, height: int = 950) -> None:
    from backend.phylogeny.viewer import load_tree_newick, newick_to_tree_json
    from backend.phylogeny.run_phylogeny import _load_backbone_metadata
    import json as _json
    _newick, _tree_err = load_tree_newick(tree_path)
    if _tree_err:
        st.warning(_tree_err)
        return
    _backbone_meta = _load_backbone_metadata()
    _tree_json_str = _json.dumps(newick_to_tree_json(
        _newick, clusters=clusters,
        upload_names=upload_names, backbone_meta=_backbone_meta,
    ))
    _html = f"""<!DOCTYPE html>
<html>
<head>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  body {{ margin: 0; background: white; font-family: Arial, sans-serif; overflow-x: hidden; }}
  .branch {{ fill: none; stroke: #444; stroke-width: 1px; }}
  .leaf-label {{ font-size: 11px; fill: #111; dominant-baseline: middle; }}
  #legend {{
    position: fixed; top: 12px; right: 12px;
    background: rgba(255,255,255,0.96); border: 1px solid #ccc;
    border-radius: 6px; padding: 8px 12px; font-size: 12px; z-index: 10;
  }}
  .legend-row {{ display: flex; align-items: center; gap: 6px; margin: 3px 0; }}
  .legend-dot {{ width: 10px; height: 10px; border-radius: 50%; flex-shrink: 0; }}
  #info-panel {{
    display: none; position: fixed; bottom: 16px; left: 16px;
    background: #fff; border: 1px solid #d1d5db; border-radius: 8px;
    padding: 12px 16px; font-size: 12px; max-width: 320px;
    box-shadow: 0 4px 12px rgba(0,0,0,0.12); z-index: 20; line-height: 1.6;
  }}
  #info-panel .ip-title {{ font-size: 14px; font-weight: 700; margin-bottom: 6px; color: #111; }}
  #info-panel .ip-row {{ display: flex; gap: 6px; color: #374151; }}
  #info-panel .ip-label {{ color: #6b7280; min-width: 70px; }}
  #info-panel .ip-resist {{ margin-top: 6px; padding: 5px 8px; background: #fef3c7;
    border-left: 3px solid #f59e0b; border-radius: 3px; color: #92400e; }}
  #info-panel .ip-close {{ float: right; cursor: pointer; color: #9ca3af;
    font-size: 16px; line-height: 1; margin-left: 8px; }}
</style>
</head>
<body>
<div id="legend"><strong>{legend_title}</strong></div>
<div id="info-panel">
  <span class="ip-close" onclick="document.getElementById('info-panel').style.display='none'">✕</span>
  <div class="ip-title" id="ip-name"></div>
  <div class="ip-row"><span class="ip-label">Origin</span><span id="ip-origin">—</span></div>
  <div class="ip-row"><span class="ip-label">MLST</span><span id="ip-mlst">—</span></div>
  <div class="ip-row"><span class="ip-label">NG-STAR</span><span id="ip-ngstar">—</span></div>
  <div class="ip-resist" id="ip-resist" style="display:none"></div>
  <div class="ip-row" id="ip-notes-row" style="margin-top:4px;color:#6b7280;font-style:italic"><span id="ip-notes"></span></div>
</div>
<svg id="tree"></svg>
<script>
const data = {_tree_json_str};
const margin = {{ top: 20, right: 460, bottom: 30, left: 50 }};
const W = (window.innerWidth || 900);
const root = d3.hierarchy(data);
const nLeaves = root.leaves().length;
const rowH = 20;
const svgH = Math.max(500, nLeaves * rowH) + margin.top + margin.bottom + 40;
const treeW = W - margin.left - margin.right;
const treeH = svgH - margin.top - margin.bottom - 40;
const svg = d3.select('#tree').attr('width', W).attr('height', svgH);
const g = svg.append('g').attr('transform', `translate(${{margin.left}},${{margin.top}})`);
svg.call(d3.zoom().scaleExtent([0.2, 8]).on('zoom', e => g.attr('transform', e.transform)));
const layout = d3.cluster().size([treeH, treeW]);
layout(root);

// Phylogram: override y positions with cumulative branch lengths
root.each(d => {{
    d.cumulLen = d.parent ? (d.parent.cumulLen || 0) + (d.data.length || 0) : 0;
}});
const maxLen = d3.max(root.descendants(), d => d.cumulLen) || 1;
root.each(d => {{ d.y = d.cumulLen / maxLen * treeW; }});

// Scale bar
const scaleBarLen = treeW * 0.1;
const scaleBarVal = (maxLen * 0.1).toFixed(4);
const sbG = g.append('g').attr('transform', `translate(0,${{treeH + 20}})`);
sbG.append('line').attr('x1', 0).attr('x2', scaleBarLen).attr('y1', 0).attr('y2', 0)
    .attr('stroke', '#666').attr('stroke-width', 1.5);
sbG.append('line').attr('x1', 0).attr('x2', 0).attr('y1', -4).attr('y2', 4)
    .attr('stroke', '#666').attr('stroke-width', 1.5);
sbG.append('line').attr('x1', scaleBarLen).attr('x2', scaleBarLen).attr('y1', -4).attr('y2', 4)
    .attr('stroke', '#666').attr('stroke-width', 1.5);
sbG.append('text').attr('x', scaleBarLen / 2).attr('y', 14)
    .attr('text-anchor', 'middle').attr('font-size', '10px').attr('fill', '#666')
    .text(scaleBarVal);

const clusterMap = {{}};
root.leaves().forEach(l => {{
    const c = l.data.cluster;
    if (c && c !== '') {{ if (!clusterMap[c]) clusterMap[c] = []; clusterMap[c].push(l); }}
}});
const shadowPalette = [
    'rgba(147,197,253,0.30)', 'rgba(249,168,212,0.30)',
    'rgba(134,239,172,0.30)', 'rgba(253,224,132,0.30)',
    'rgba(196,181,253,0.30)', 'rgba(103,232,249,0.30)',
    'rgba(254,202,202,0.30)', 'rgba(167,243,208,0.30)',
];
let ci = 0;
Object.entries(clusterMap).filter(([, ls]) => ls.length >= 2).sort(([a],[b]) => +a - +b)
    .forEach(([c, ls]) => {{
        const xs = ls.map(l => l.x);
        const yMin = Math.min(...xs) - rowH * 0.6;
        const yMax = Math.max(...xs) + rowH * 0.6;
        const col = shadowPalette[ci++ % shadowPalette.length];
        g.append('rect').attr('x', 0).attr('y', yMin)
            .attr('width', treeW + 10).attr('height', yMax - yMin)
            .attr('rx', 8).attr('fill', col);
        const lx = treeW + 240, ly = (yMin + yMax) / 2, lh = 22;
        g.append('rect').attr('x', lx - 8).attr('y', ly - lh/2)
            .attr('width', 110).attr('height', lh).attr('rx', 5)
            .attr('fill', col.replace('0.30','0.70')).attr('stroke', col.replace('0.30','0.90'))
            .attr('stroke-width', 1);
        g.append('text').attr('x', lx).attr('y', ly).attr('dominant-baseline', 'middle')
            .attr('font-size', '12px').attr('font-weight', '700').attr('fill', '#1f2937')
            .text('Cluster ' + c);
    }});
g.selectAll('.branch').data(root.links()).join('path').attr('class', 'branch')
    .attr('d', d => `M${{d.source.y}},${{d.source.x}}L${{d.source.y}},${{d.target.x}}L${{d.target.y}},${{d.target.x}}`);
g.selectAll('.inode').data(root.descendants().filter(d => d.children)).join('circle')
    .attr('cx', d => d.y).attr('cy', d => d.x).attr('r', 2).attr('fill', '#aaa');
const leaves = g.selectAll('.leaf').data(root.leaves()).join('g')
    .attr('class', 'leaf').attr('transform', d => `translate(${{d.y}},${{d.x}})`);
leaves.append('circle').attr('r', 5).attr('fill', d => d.data.color)
    .attr('stroke', d => d.data.type === 'leaf_upload' ? '#333' : 'none').attr('stroke-width', 1.5);
leaves.append('text').attr('class', 'leaf-label').attr('x', 9).text(d => d.data.name);

// Click handler — show metadata panel for any leaf
leaves.style('cursor', 'pointer').on('click', function(event, d) {{
    event.stopPropagation();
    const info = d.data.backbone_info || null;
    const panel = document.getElementById('info-panel');
    document.getElementById('ip-name').textContent =
        info ? (info.label || d.data.name) : d.data.name;
    const origin = info
        ? [info.country, info.year].filter(Boolean).join(' · ')
        : (d.data.type === 'leaf_upload' ? 'Uploaded sample' : '—');
    document.getElementById('ip-origin').textContent = origin || '—';
    document.getElementById('ip-mlst').textContent   = (info && info.mlst_st)   || '—';
    document.getElementById('ip-ngstar').textContent = (info && info.ngstar_st) || '—';
    const resistEl = document.getElementById('ip-resist');
    if (info && info.resistance_profile) {{
        resistEl.textContent = info.resistance_profile;
        resistEl.style.display = 'block';
    }} else {{
        resistEl.style.display = 'none';
    }}
    const notesEl  = document.getElementById('ip-notes');
    const notesRow = document.getElementById('ip-notes-row');
    if (info && info.notes) {{
        notesEl.textContent  = info.notes;
        notesRow.style.display = 'flex';
    }} else {{
        notesRow.style.display = 'none';
    }}
    panel.style.display = 'block';
}});
// Click on background closes panel
document.getElementById('tree').addEventListener('click', () => {{
    document.getElementById('info-panel').style.display = 'none';
}});
const legend = document.getElementById('legend');
const typesSeen = {{}};
root.leaves().forEach(n => {{
    const key = n.data.type;
    const label = key === 'leaf_upload' ? 'Uploaded' : key === 'outgroup' ? 'Outgroup' : 'Reference';
    if (!typesSeen[key]) {{
        typesSeen[key] = true;
        const row = document.createElement('div');
        row.className = 'legend-row';
        row.innerHTML = `<span class="legend-dot" style="background:${{n.data.color}}"></span>${{label}}`;
        legend.appendChild(row);
    }}
}});
</script>
</body>
</html>"""
    st.components.v1.html(_html, height=height, scrolling=True)


def render() -> None:
    _hdr_col, _report_col, _share_col = st.columns([3, 1, 1])
    with _hdr_col:
        st.markdown("""
<div style="padding:0.6rem 0 0.2rem 0;">
  <div style="font-size:0.7rem;font-weight:700;text-transform:uppercase;letter-spacing:0.1em;
              color:#64748b;margin-bottom:0.25rem;">Analysis Report</div>
  <div style="font-size:1.6rem;font-weight:800;color:#0f172a;line-height:1.2;">Full Pipeline Results</div>
</div>
""", unsafe_allow_html=True)

    _fpr = st.session_state.get("full_pipeline_results")
    if not _fpr:
        st.info("No pipeline results available yet. Run the Full Pipeline from the Homepage.")
        st.stop()

    _active_proj = st.session_state.get("active_project", "")
    if _active_proj and db:
        with _report_col:
            st.markdown("<div style='padding-top:0.55rem'></div>", unsafe_allow_html=True)
            try:
                from backend.report import generate_html_report
                _db_for_report = _load_project_cached(_active_proj)
                _report_bytes  = generate_html_report(
                    _active_proj,
                    _db_for_report["samples"],
                    _db_for_report.get("phylogeny"),
                ).encode("utf-8")
                st.download_button(
                    "⬇ Report (HTML)",
                    data=_report_bytes,
                    file_name=f"{_active_proj}_report.html",
                    mime="text/html",
                    use_container_width=True,
                    key="report_dl_btn",
                )
            except Exception as _rep_err:
                st.caption(f"Report error: {_rep_err}")

    if _active_proj:
        with _share_col:
            st.markdown("<div style='padding-top:0.55rem'></div>", unsafe_allow_html=True)
            if st.button("Share results", use_container_width=True):
                st.session_state["_show_share"] = not st.session_state.get("_show_share", False)

        if st.session_state.get("_show_share"):
            _safe_proj = _active_proj.replace("'", "\\'")
            st.components.v1.html(f"""
<div style="display:flex;gap:8px;align-items:center;padding:2px 0 10px 0;">
  <input id="ng-share-url" type="text" value="Generating link…"
    style="flex:1;padding:8px 14px;border:1.5px solid #e5e7eb;border-radius:10px;
           font-size:0.82rem;font-family:'Courier New',monospace;background:#f9fafb;
           color:#111;outline:none;"
    readonly onclick="this.select()">
  <button id="ng-copy-btn"
    onclick="navigator.clipboard.writeText(document.getElementById('ng-share-url').value)
             .then(()=>{{var b=document.getElementById('ng-copy-btn');b.textContent='✓ Copied!';
               b.style.background='#16a34a';setTimeout(()=>{{b.textContent='Copy link';b.style.background='#1d4ed8';}},2000)}})
             .catch(()=>{{document.getElementById('ng-share-url').select();document.execCommand('copy');}})"
    style="padding:8px 18px;background:#1d4ed8;color:white;border:none;border-radius:10px;
           cursor:pointer;font-size:0.82rem;font-weight:600;white-space:nowrap;transition:background .2s;">
    Copy link
  </button>
</div>
<div style="font-size:0.7rem;color:#6b7280;margin-top:-4px;padding-bottom:4px;">
  Anyone with access to this server can open this link and view these results directly.
</div>
<script>
(function() {{
  try {{
    var url = new URL(window.location.href);
    url.searchParams.set('project', '{_safe_proj}');
    url.searchParams.set('page', 'results');
    document.getElementById('ng-share-url').value = url.toString();
  }} catch(e) {{
    document.getElementById('ng-share-url').value = window.location.origin + '/?project={_safe_proj}&page=results';
  }}
}})();
</script>
""", height=80)

   
    if db and _active_proj:
        try:
            from backend.alerts import detect_alerts as _detect_alerts
            _db_snap = _load_project_cached(_active_proj)
            _new_alerts = _detect_alerts(_db_snap["samples"])
            db.save_alerts(_active_proj, _new_alerts)

            _all_alerts = db.load_alerts(_active_proj)
            _unread     = [a for a in _all_alerts if not a["acknowledged"]]

            if _all_alerts:
                _sev_order = {"critical": 0, "error": 1, "warning": 2}
                _unread_sorted = sorted(_unread, key=lambda a: _sev_order.get(a["severity"], 9))

                _sev_icons = {"critical": "🚨", "error": "🔴", "warning": "🟡"}
                _sev_colors = {"critical": "#701a75", "error": "#991b1b", "warning": "#92400e"}

                from collections import Counter as _Counter
                _counts = _Counter(a["severity"] for a in _unread)
                _parts = [
                    f"{_sev_icons[s]} {_counts[s]} {s}"
                    for s in ("critical", "error", "warning") if _counts.get(s)
                ]
                _exp_label = " · ".join(_parts) if _parts else "✅ No new alerts"

                with st.expander(_exp_label, expanded=False):
                    if _unread_sorted:
                        for _al in _unread_sorted:
                            _icon = _sev_icons.get(_al["severity"], "ℹ️")
                            _fg   = _sev_colors.get(_al["severity"], "#0369a1")
                            st.markdown(
                                f'<div style="padding:4px 0;font-size:13px;">'
                                f'<span style="font-weight:600;color:{_fg}">{_icon} {_al["message"]}</span>'
                                f'</div>',
                                unsafe_allow_html=True,
                            )
                        st.markdown("")
                        if st.button("Mark all as read", key="ack_all_alerts"):
                            db.acknowledge_alerts(_active_proj)
                            st.rerun()
                    else:
                        st.success("All alerts acknowledged.")

                    if len(_all_alerts) > len(_unread):
                        with st.expander("Previously acknowledged alerts"):
                            for _al in [a for a in _all_alerts if a["acknowledged"]]:
                                st.caption(f"~~{_al['message']}~~ _{_al['created_at'][:16]}_")
        except Exception:
            pass

    _fpr_sids = sorted(set(
        list(_fpr.get("qc", {}).keys())
        + list(_fpr.get("assembly", {}).keys())
        + list(_fpr.get("amr", {}).keys())
    ))

    _fpr_meta: dict = {}
    if _active_proj and db:
        try:
            _fpr_meta = db.load_metadata(_active_proj)
        except Exception:
            pass

    #  Genomic overview 
    st.markdown("""
<div style="margin:1.2rem 0 0.8rem 0;padding-bottom:0.55rem;border-bottom:2px solid #e2e8f0;">
  <div style="font-size:1.05rem;font-weight:700;color:#0f172a;">Genomic Overview</div>
</div>
""", unsafe_allow_html=True)
    _ov_rows = []
    for _sid in _fpr_sids:
        _qr = _fpr.get("qc", {}).get(_sid, {})
        _ar = _fpr.get("assembly", {}).get(_sid, {})
        _mr = _fpr.get("amr", {}).get(_sid, {})

        _m  = _qr.get("metrics", {}) if not _qr.get("error") else {}
        _k2 = _qr.get("kraken2", {}) if not _qr.get("error") else {}
        _sp = _k2.get("species_status", "—")
        _sp_icon = {"confirmed": "✅", "likely": "⚠️", "contaminated": "❌"}.get(_sp, "—")

        _ctg_path = _ar.get("contigs_path") if not _ar.get("error") else None
        _cst = parse_contigs_stats(_ctg_path) if _ctg_path else None
        _asm_qc_icon = {"pass": "✅", "warn": "⚠️", "fail": "❌"}.get(
            _cst.qc_status if _cst is not None else None, "—"
        )

        _prob   = _mr.failure_probability if isinstance(_mr, AMRResult) else None
        _cdc    = _mr.cdc_phenotypes if isinstance(_mr, AMRResult) else []
        _pheno  = "; ".join(c.replace("_", " ").capitalize() for c in _cdc if c != "wildtype") or "Wildtype"
        _mlst   = _mr.mlst if isinstance(_mr, AMRResult) else {}
        _ngstar = _mr.ngstar if isinstance(_mr, AMRResult) else {}
        _st_str     = f"ST-{_mlst['st']}"   if not _mlst.get("error")   and _mlst.get("st")  else "—"
        _ngstar_str = f"ST-{_ngstar['ST']}" if not _ngstar.get("error") and _ngstar.get("ST") else "—"
        _res_cat     = _mr.resistance_category if isinstance(_mr, AMRResult) else None
        _res_cat_str = _res_cat.replace("_", " ").capitalize() if _res_cat else "—"

        _cov   = _m.get("mean_coverage")
        _b10x  = _m.get("pct_breadth_10x")
        _cov_str  = f"{_cov:.1f}×" if _cov is not None else "—"
        _b10x_str = f"{_b10x:.1f}%" if _b10x is not None else "—"

        _smeta = _fpr_meta.get(_sid, {})
        _ov_rows.append({
            "Sample":               _sid,
            "Date":                 _smeta.get("collection_date") or "—",
            "Site":                 _smeta.get("anatomical_site") or "—",
            "Location":             ", ".join(filter(None, [_smeta.get("city"), _smeta.get("country")])) or "—",
            "Reads":                f"{_m['total_reads']:,}" if _m.get("total_reads") else "—",
            "Coverage":             _cov_str,
            "≥10×":                 _b10x_str,
            "Species":              f"{_sp_icon} {_sp.replace('_',' ').capitalize()}" if _sp != "—" else "—",
            "Assembly QC":          _asm_qc_icon,
            "Assembly":             f"{_cst.total_len/1e6:.2f} Mb" if _cst is not None else "—",
            "N50":                  f"{_cst.n50/1e3:.1f} kb" if _cst is not None else "—",
            "Core genes %":         (
                f"{_cst.completeness:.1f}% ({_cst.core_genes_found}/{_cst.core_genes_total})"
                if _cst is not None and _cst.core_genes_found is not None
                else (f"{_cst.completeness:.1f}%" if _cst is not None else "—")
            ),
            "ST (MLST)":            _st_str,
            "NG-STAR":              _ngstar_str,
            "Resistance category":  _res_cat_str,
            "AMR score":            f"{_prob*100:.1f}%" if _prob is not None else "—",
            "CDC phenotype":        _pheno,
        })
    _ov_event = st.dataframe(
        pd.DataFrame(_ov_rows),
        use_container_width=True,
        hide_index=True,
        selection_mode="single-row",
        on_select="rerun",
        key="fp_ov_tbl",
    )

    _ov_selected = (
        _ov_event.selection.rows
        if _ov_event and hasattr(_ov_event, "selection")
        else []
    )

    if _ov_selected:
        _sel_sid = _fpr_sids[_ov_selected[0]]
        st.markdown(f"#### {_sel_sid} — full report")

        _s_qr = _fpr.get("qc",       {}).get(_sel_sid, {})
        _s_ar = _fpr.get("assembly",  {}).get(_sel_sid, {})
        _s_mr = _fpr.get("amr",       {}).get(_sel_sid, {})

        _sel_meta = _fpr_meta.get(_sel_sid, {})
        if any(_sel_meta.values()):
            _mc1, _mc2, _mc3, _mc4 = st.columns(4)
            _mc1.metric("Collection date", _sel_meta.get("collection_date") or "—")
            _mc2.metric("Site", _sel_meta.get("anatomical_site") or "—")
            _mc3.metric("Location", _sel_meta.get("city") or _sel_meta.get("country") or "—")
            _mc4.metric("Health unit", _sel_meta.get("health_unit") or "—")

        _rc1, _rc2, _rc3 = st.columns(3)

        with _rc1:
            st.markdown("**QC**")
            _sm = _s_qr.get("metrics", {})
            if _sm:
                st.write(f"Reads: {_sm.get('total_reads', '—'):,}" if _sm.get('total_reads') else "Reads: —")
                st.write(f"GC: {_sm.get('gc', '—')}%" if _sm.get('gc') else "GC: —")
                _cov = _sm.get('mean_coverage')
                st.write(f"Coverage: {_cov:.1f}×" if _cov else "Coverage: —")
            else:
                st.write(_s_qr.get("error", "—"))

        with _rc2:
            st.markdown("**Assembly**")
            _s_cst = parse_contigs_stats(_s_ar.get("contigs_path")) if _s_ar.get("contigs_path") else None
            if _s_cst is not None:
                st.write(f"Length: {_s_cst.total_len/1e6:.2f} Mb")
                st.write(f"N50: {_s_cst.n50/1e3:.1f} kb")
                st.write(f"Contigs: {_s_cst.n_contigs}")
                st.write(f"GC: {_s_cst.gc_pct:.1f}%")
                st.write(f"Assembly QC: {_s_cst.qc_status}")
            else:
                st.write(_s_ar.get("error", "—"))

        with _rc3:
            st.markdown("**Download**")
            _s_ctg = _s_ar.get("contigs_path")
            if _s_ctg and Path(_s_ctg).exists():
                st.download_button(
                    "⬇ Contigs FASTA",
                    data=Path(_s_ctg).read_bytes(),
                    file_name=f"{_sel_sid}_contigs.fasta",
                    mime="application/octet-stream",
                    key=f"fp_det_dl_{_sel_sid}",
                )
            _s_mqc = Path(_s_qr.get("multiqc", "") or "")
            if _s_mqc.is_file():
                st.download_button(
                    "⬇ MultiQC report",
                    data=_s_mqc.read_bytes(),
                    file_name=f"{_sel_sid}_multiqc.html",
                    mime="text/html",
                    key=f"fp_mqc_det_{_sel_sid}",
                )

        if isinstance(_s_mr, AMRResult):
            _det_mlst   = _s_mr.mlst
            _det_ngstar = _s_mr.ngstar
            _det_st     = f"ST-{_det_mlst['st']}"   if not _det_mlst.get("error")   and _det_mlst.get("st")  else "—"
            _det_ng     = f"ST-{_det_ngstar['ST']}"  if not _det_ngstar.get("error") and _det_ngstar.get("ST") else "—"
            _tc1, _tc2 = st.columns(2)
            _tc1.metric("MLST", _det_st)
            _tc2.metric("NG-STAR", _det_ng)
            st.markdown("**AMR Profile**")
            render_amr_interpretation(_s_mr, sample_label=_sel_sid, contigs_path=_s_ctg)
        elif isinstance(_s_mr, dict) and _s_mr.get("error"):
            st.error(f"AMR error: {_s_mr['error']}")
    else:
        st.caption("Click a row to view the full report for that sample.")

    st.divider()

    #  AMR profiling
    st.markdown("""
<div style="margin:1.2rem 0 0.8rem 0;padding-bottom:0.55rem;border-bottom:2px solid #e2e8f0;">
  <div style="font-size:1.05rem;font-weight:700;color:#0f172a;">AMR Profiling</div>
</div>
""", unsafe_allow_html=True)

    _amr_data  = _fpr.get("amr", {})
    _amr_valid = {sid: r for sid, r in _amr_data.items() if isinstance(r, AMRResult)}

    if _amr_valid:
        render_amr_mutation_matrix(_amr_valid)
        _render_essential_gene_markers(_amr_valid)

        st.markdown("#### Clinical Interpretation")
        for _sid, _res in _amr_valid.items():
            _prob  = _res.failure_probability
            _cdc   = _res.cdc_phenotypes
            _ther  = _res.therapy
            _mlst  = _res.mlst
            _mos   = _res.mosaic_pena
            _is_wt = _cdc == ["wildtype"]

            if _is_wt:
                _sbg, _sc, _sb = "#f0fdf4", "#16a34a", "#bbf7d0"
                _slabel, _sicon = "Wildtype — no resistance detected", "✅"
            else:
                _sbg, _sc = _resistance_color(_prob)
                _sb = _sbg
                _n_r = len([p for p in _cdc if p != "wildtype"])
                _slabel = f"{_n_r} resistance mechanism(s) detected"
                _sicon = "🔴" if _prob >= 0.4 else "🟡"

            _st_s = f"ST-{_mlst['st']}" if not _mlst.get("error") and _mlst.get("st") else "—"
            _mos_s = "Yes" if _mos.get("mosaic_suspected") else ("No" if _mos.get("identity") else "—")

            with st.expander(
                f"**{_sid}** — {_sicon} {_slabel} · Score {_prob*100:.1f}%",
                expanded=len(_amr_valid) <= 3,
            ):
                st.caption(f"MLST: {_st_s}  ·  Mosaic penA: {_mos_s}  ·  European 2020 (IUSTI) guidelines")
                _avoid = _ther.get("avoid", [])
                _rec   = _ther.get("recommend", [])
                _avoid_low_fp = [a.lower() for a in _avoid]
                _rec_str_fp   = " ".join(_rec).lower()
                _fp_drug_cards = ""
                for _flabel, _fkey in [
                    ("Ceftriaxone", "ceftriaxone"), ("Azithromycin", "azithromycin"),
                    ("Gentamicin", "gentamicin"), ("Ciprofloxacin", "ciprofloxacin"),
                    ("Penicillin", "penicillin"), ("Tetracycline", "tetracycline"),
                    ("Doxycycline", "doxycycline"),
                ]:
                    if _fkey in _avoid_low_fp:
                        _fbg, _fbrd, _ftxt, _ficon, _fsub = "#fef2f2", "#dc2626", "#991b1b", "✕", "Avoid"
                    elif _fkey in _rec_str_fp:
                        _fbg, _fbrd, _ftxt, _ficon, _fsub = "#f0fdf4", "#16a34a", "#15803d", "✓", "Recommended"
                    else:
                        _fbg, _fbrd, _ftxt, _ficon, _fsub = "#f9fafb", "#d1d5db", "#6b7280", "—", "Not indicated"
                    _fp_drug_cards += (
                        f'<div style="background:{_fbg};border:1.5px solid {_fbrd};border-radius:10px;'
                        f'padding:0.7rem 0.8rem;text-align:center;min-width:90px;flex:1;">'
                        f'<div style="font-size:1.1rem;font-weight:700;color:{_ftxt};">{_ficon}</div>'
                        f'<div style="font-size:0.78rem;font-weight:700;color:{_ftxt};margin:0.15rem 0;">{_flabel}</div>'
                        f'<div style="font-size:0.65rem;color:{_ftxt};opacity:0.8;">{_fsub}</div>'
                        f'</div>'
                    )
                st.markdown(
                    f'<div style="display:flex;gap:0.4rem;flex-wrap:wrap;margin-bottom:0.8rem;">{_fp_drug_cards}</div>',
                    unsafe_allow_html=True,
                )
                _fp_ctg = (_fpr.get("assembly") or {}).get(_sid, {}).get("contigs_path")
                _render_genome_overview(_fp_ctg or "", _res, _sid)

        for _sid, _err in {s: r["error"] for s, r in _amr_data.items() if "error" in r}.items():
            st.error(f"**{_sid}**: {_err}")
        if _active_proj and db:
            _fp_cohort_amr = {
                sid: d["amr"]
                for sid, d in _load_project_cached(_active_proj)["samples"].items()
                if "amr" in d
            }
            _render_cohort_amr_profile(_fp_cohort_amr, project_name=_active_proj)
    else:
        st.info("AMR profiling not yet completed.")

    st.divider()

    st.markdown("""
<div style="margin:1.2rem 0 0.8rem 0;padding-bottom:0.55rem;border-bottom:2px solid #e2e8f0;">
  <div style="font-size:1.05rem;font-weight:700;color:#0f172a;">Fast NJ Phylogeny (SKA2)</div>
</div>
""", unsafe_allow_html=True)
    _phy = _fpr.get("phylogeny", {})

    if _phy.get("error"):
        st.warning(f"Phylogeny: {_phy['error']}")
    else:
        _tree_path = st.session_state.get("last_tree_path") or _phy.get("tree_path")
        _clusters  = st.session_state.get("last_clusters") or _phy.get("clusters", {})
        _cr        = st.session_state.get("last_cluster_report") or _phy.get("cluster_report", {})

        if _tree_path and Path(_tree_path).exists():
            if _clusters:
                from collections import Counter as _Counter
                _ctr = _Counter(_clusters.values())
                _n_multi = sum(1 for v in _ctr.values() if v >= 2)
                st.caption(f"NJ clusters (SKA2): {len(_ctr)} clusters · {_n_multi} with ≥2 members")
            _upload_names = st.session_state.get("last_upload_names", [])
            try:
                _render_d3_tree(_tree_path, _clusters, _upload_names, "NJ clusters (SKA2)")
            except Exception as _te:
                st.warning(f"Could not render tree: {_te}")
        else:
            st.info("Tree not yet available.")

        if _cr:
            _up_rows = []
            for _n, _i in sorted(_cr.items()):
                if not _i.get("is_uploaded"):
                    continue
                _nn      = _i.get("nearest_neighbour", "—")
                _nn_info = _i.get("nn_backbone_info", {})
                _nn_label = _nn_info.get("label", _nn) if _nn_info else _nn
                _nn_origin = ""
                if _nn_info:
                    _parts = [_nn_info.get("country"), str(_nn_info.get("year", ""))]
                    _nn_origin = " · ".join(p for p in _parts if p and p != "None")
                _up_rows.append({
                    "Sample":            _n,
                    "Cluster":           str(_i.get("cluster", "—")),
                    "Nearest neighbour": _nn_label,
                    "Origin":            _nn_origin or "—",
                    "Distance (SNPs)":   f"{_i.get('nn_distance', 0):.0f}",
                    "Distance (%)":      f"{_i.get('nn_distance_pct', 0.0):.4f}%",
                    "NN resistance":     _nn_info.get("resistance_profile", "—") if _nn_info else "—",
                })
            if _up_rows:
                st.markdown("**Uploaded samples — cluster assignment:**")
                st.dataframe(pd.DataFrame(_up_rows), use_container_width=True, hide_index=True)

            _nn_cards = [
                (_n, _i) for _n, _i in sorted(_cr.items())
                if _i.get("is_uploaded") and _i.get("nearest_neighbour")
            ]
            if _nn_cards:
                st.markdown("**Closest reference genomes:**")
                _card_html_parts = []
                for _n, _i in _nn_cards:
                    _nn      = _i.get("nearest_neighbour", "—")
                    _nn_info = _i.get("nn_backbone_info", {})
                    _snps    = _i.get("nn_distance", 0)
                    _pct     = _i.get("nn_distance_pct", 0.0)

                    if _nn_info:
                        _nn_label  = _nn_info.get("label", _nn)
                        _nn_origin = " · ".join(
                            str(x) for x in [_nn_info.get("country"), _nn_info.get("year")]
                            if x and str(x) != "None"
                        )
                        _nn_resist  = _nn_info.get("resistance_profile", "")
                        _nn_classes = _nn_info.get("resistance_classes", [])
                        _nn_mlst    = _nn_info.get("mlst_st") or ""

                        if not _nn_classes:
                            _card_bg, _card_brd, _resist_bg, _resist_col = (
                                "#f0fdf4", "#86efac", "#dcfce7", "#15803d")
                            _resist_icon = "✓"
                        elif any(c in _nn_classes for c in ("ceftriaxone", "azithromycin")):
                            _card_bg, _card_brd, _resist_bg, _resist_col = (
                                "#fef2f2", "#fca5a5", "#fee2e2", "#dc2626")
                            _resist_icon = "⚠"
                        else:
                            _card_bg, _card_brd, _resist_bg, _resist_col = (
                                "#fffbeb", "#fcd34d", "#fef3c7", "#b45309")
                            _resist_icon = "⚠"

                        _origin_str = f"<span style='color:#6b7280'>{_nn_origin}</span>" if _nn_origin else ""
                        _mlst_str   = f"<span style='color:#6b7280;margin-left:8px;'>{_nn_mlst}</span>" if _nn_mlst else ""
                        _resist_str = (
                            f"<div style='margin-top:5px;padding:3px 8px;background:{_resist_bg};"
                            f"border-left:3px solid {_resist_col};border-radius:3px;"
                            f"font-size:0.78rem;color:{_resist_col};'>"
                            f"{_resist_icon} {_nn_resist}</div>"
                        ) if _nn_resist else ""
                    else:
                        _nn_label   = _nn
                        _origin_str = "<span style='color:#6b7280'>uploaded sample</span>"
                        _mlst_str   = ""
                        _resist_str = ""
                        _card_bg, _card_brd = "#f9fafb", "#d1d5db"

                    _card_html_parts.append(f"""
<div style="border:1px solid {_card_brd};background:{_card_bg};border-radius:8px;
            padding:8px 12px;margin-bottom:6px;font-size:0.85rem;">
  <div style="display:flex;justify-content:space-between;align-items:baseline;">
    <span><strong>{_n}</strong>
      <span style="color:#9ca3af;margin:0 6px;">→</span>
      <strong>{_nn_label}</strong>
      &nbsp;{_origin_str}{_mlst_str}
    </span>
    <span style="color:#374151;font-size:0.78rem;white-space:nowrap;margin-left:12px;">
      {int(_snps)} SNPs &nbsp;·&nbsp; {_pct:.4f}%
    </span>
  </div>
  {_resist_str}
</div>""")

                st.markdown(
                    "<div style='margin-top:4px'>" + "".join(_card_html_parts) + "</div>",
                    unsafe_allow_html=True,
                )



