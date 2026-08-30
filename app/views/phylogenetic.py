from __future__ import annotations
import json
import time
from pathlib import Path

import pandas as pd
import streamlit as st

from views.shared import *
from views.shared import _inline_metadata_widget, _render_cgmlst_clustering, _render_amr_tree_highlight


def render() -> None:

    st.html("""
<div style="padding:0.6rem 0 1rem 0;">
  <div style="display:inline-flex;align-items:center;background:rgba(0,0,0,0.04);border:1px solid rgba(0,0,0,0.12);
  border-radius:20px;padding:0.22rem 0.85rem;font-size:0.7rem;font-weight:700;color:#64748b;
  letter-spacing:0.12em;text-transform:uppercase;margin-bottom:0.6rem;">Module 3</div>
  <div style="font-size:1.6rem;font-weight:800;color:#0f172a;">3. Phylogenetic Analysis</div>
</div>
""")

    _phy_active_proj = st.session_state.get("active_project")
    if _phy_active_proj and db:
        if st.button(" Delete all results", key="del_phy_btn"):
            from backend.db import delete_phylogeny as _delete_phylogeny
            _delete_phylogeny(_phy_active_proj)
            for _k in ("phy_run_id", "phy_tree_path", "phy_tree_method",
                       "phy_cluster_report", "phy_matrix_path", "phy_upload_names",
                       "phy_cgmlst_result",
                       "last_tree_path", "last_clusters", "last_cluster_report",
                       "last_matrix_path", "last_upload_names", "last_cgmlst_result"):
                st.session_state.pop(_k, None)
            _load_project_cached.clear()
            st.rerun()

    from backend.phylogeny.viewer import load_tree_newick, newick_to_tree_json

    from backend.phylogeny.cgmlst import is_schema_ready, download_schema
    from backend.phylogeny.run_pyngost import is_pyngost_ready, download_pyngost_db
    _schema_ready  = is_schema_ready()
    _pyngost_ready = is_pyngost_ready()
    if not _schema_ready or not _pyngost_ready:
        with st.expander("Reference databases missing", expanded=True):
            if not _schema_ready:
                st.warning(
                    "cgMLST schema not found "
                )
                if st.button("Download cgMLST schema", key="dl_cgmlst_schema"):
                    with st.status("Downloading cgMLST schema from PubMLST…", expanded=True) as _s:
                        try:
                            download_schema()
                            _s.update(label="cgMLST schema ready.", state="complete")
                            st.rerun()
                        except Exception as e:
                            _s.update(label="Download failed", state="error")
                            st.error(str(e))
            if not _pyngost_ready:
                st.warning(
                    "pyngoST allele database not found."
                )
                if st.button("Download pyngoST database", key="dl_pyngost_db"):
                    with st.status("Downloading MLST/NG-STAR databases…", expanded=True) as _s:
                        try:
                            download_pyngost_db()
                            _s.update(label="pyngoST database ready.", state="complete")
                            st.rerun()
                        except Exception as e:
                            _s.update(label="Download failed", state="error")
                            st.error(str(e))


    st.markdown("""
    This module builds WGS phylogenetic trees for *Neisseria gonorrhoeae* isolates,
    placed in context with a curated backbone panel of international reference genomes and
    Portuguese clinical isolates. Accepts assembled genomes (FASTA). Raw FASTQ reads
    can be uploaded in Module 1 (Quality Control & Assembly) or via the Full Pipeline in HomePage.

    1. **SKA2** — computes pairwise SNP (allelic) distances between all genomes
    2. **RapidNJ** — builds the Neighbour-Joining tree from the distance matrix
    3. **cgMLST clustering** (chewBBACA) — genogroup clusters assigned at a fixed threshold
       of ≤400 allele differences (AD)
    """)

    uploaded_genomes = st.file_uploader(
        "Upload assembled genomes (FASTA).",
        accept_multiple_files=True,
        type=["fasta", "fa", "fna"],
        key="phy_upload_fasta",
    ) or []

    _phy_all_sids = [Path(f.name).stem for f in uploaded_genomes]
    if _phy_all_sids:
        _phy_proj = st.session_state.get("active_project")
        _inline_metadata_widget(_phy_all_sids, _phy_proj, "phy")

    _phy_has_input = bool(uploaded_genomes)
    _phy_busy = bool(st.session_state.get("_phy_pending_job"))
    if _phy_has_input and st.button("Run Fast NJ", disabled=_phy_busy):
        from backend.paths import RUNS_DIR as _RUNS_DIR
        import uuid as _uuid
        _stage = Path(_RUNS_DIR) / f"stage_{_uuid.uuid4().hex}"
        _stage.mkdir(parents=True, exist_ok=True)

        _samples_payload = []
        _amr_targets: list[tuple[str, str]] = []
        for _f in uploaded_genomes:
            _dest = _stage / _f.name
            _dest.write_bytes(_f.read())
            _samples_payload.append({"kind": "fasta", "path": str(_dest)})
            _amr_targets.append((Path(_f.name).stem, str(_dest)))

        _jid = db.submit_job(
            "_phylogeny", "_phy_nj", "phylogeny_nj", {"samples": _samples_payload}
        )
        st.session_state["_phy_pending_job"] = _jid

        _phy_amr_proj = st.session_state.get("active_project")
        if _phy_amr_proj and db:
            st.session_state["_phy_amr_pending_jobs"] = [
                db.submit_job(_phy_amr_proj, sid, "amr", {"contigs": cpath})
                for sid, cpath in _amr_targets
            ]
            st.session_state["_phy_amr_pending_project"] = _phy_amr_proj
        st.rerun()

    @st.fragment(run_every="2s")
    def _phy_poll_status():
        _phy_pending_job = st.session_state.get("_phy_pending_job")
        if not (_phy_pending_job and db):
            return
        _phy_all_jobs = {j["id"]: j for j in db.get_project_jobs("_phylogeny")}
        _phy_job      = _phy_all_jobs.get(_phy_pending_job)
        if not _phy_job:
            return
        _phy_jstatus = _phy_job["status"]
        if _phy_jstatus in ("queued", "running"):
            if _phy_jstatus == "queued":
                st.info("Queued…")
            else:
                if not st.session_state.get("_phy_running_since"):
                    st.session_state["_phy_running_since"] = time.time()
                _phy_elapsed = time.time() - st.session_state["_phy_running_since"]
                if _phy_elapsed < 15:
                    _phy_stage_msg = "Building distance matrix (SKA2)"
                elif _phy_elapsed < 30:
                    _phy_stage_msg = "Building the tree (RapidNJ) and assigning clusters"
                else:
                    _phy_stage_msg = "Running cgMLST clustering"
                with st.status("Running phylogenetic analysis…", state="running", expanded=True):
                    st.write(_phy_stage_msg)
                    st.write(
                        "cgMLST step is slow and can take several minutes, but this page updates "
                        "automatically, no need to resubmit."
                    )
        elif _phy_jstatus == "done" and _phy_job.get("result"):
            _phy_res  = json.loads(_phy_job["result"])
            st.session_state["phy_run_id"]         = _phy_res["run_id"]
            st.session_state["phy_tree_path"]      = _phy_res["tree_path"]
            st.session_state["phy_tree_method"]    = "Fast NJ"
            st.session_state["phy_cluster_report"] = _phy_res.get("cluster_report", {})
            st.session_state["phy_matrix_path"]    = _phy_res.get("matrix_path")
            st.session_state["phy_upload_names"]   = _phy_res.get("upload_names", [])
            st.session_state["phy_cgmlst_result"]  = _phy_res.get("cgmlst_result")
            st.session_state.pop("_phy_pending_job", None)
            st.session_state.pop("_phy_running_since", None)
            st.rerun()
        elif _phy_jstatus == "error":
            st.error(f"Phylogeny error: {_phy_job.get('error_msg', 'Unknown error')}")
            st.session_state.pop("_phy_pending_job", None)
            st.session_state.pop("_phy_running_since", None)

    if st.session_state.get("_phy_pending_job") and db:
        _phy_poll_status()

    @st.fragment(run_every="2s")
    def _phy_amr_poll_status():
        _amr_ids = st.session_state.get("_phy_amr_pending_jobs", [])
        if not (_amr_ids and db):
            return
        _amr_proj = st.session_state.get("_phy_amr_pending_project", "")
        _all_amr_jobs = {j["id"]: j for j in db.get_project_jobs(_amr_proj)}
        _rel_amr = [_all_amr_jobs[jid] for jid in _amr_ids if jid in _all_amr_jobs]
        _n_done = sum(1 for j in _rel_amr if j["status"] in ("done", "error"))
        if _n_done < len(_rel_amr):
            st.caption(f"Computing AMR for resistance highlighting: {_n_done}/{len(_rel_amr)} done…")
        else:
            st.session_state["_phy_amr_pending_jobs"] = []
            _load_project_cached.clear()
            st.rerun()

    if st.session_state.get("_phy_amr_pending_jobs") and db:
        _phy_amr_poll_status()

    st.subheader("Tree visualization")

    tree_path = st.session_state.get("phy_tree_path", None)

    if not tree_path:
        st.info("No tree generated yet, upload genomes above and click Run.")

    else:
        newick, err = load_tree_newick(tree_path)

        if err:
            st.warning(err)

        else:
            _tree_method = st.session_state.get("phy_tree_method", "Fast NJ")
            st.success(f"Tree loaded, generated with: **{_tree_method}**")

            _cgmlst_for_tree = st.session_state.get("phy_cgmlst_result") or {}
            clusters = _cgmlst_for_tree.get("genogroup_clusters") or {}
            cluster_source = "cgMLST genogroup (≤400 AD)"
            upload_names = st.session_state.get("phy_upload_names", [])

            if clusters:
                from collections import Counter
                ctr = Counter(clusters.values())
                n_multi = sum(1 for v in ctr.values() if v >= 2)
                st.caption(f"{len(ctr)} clusters detected ({cluster_source}), {n_multi} with ≥2 members")
            else:
                st.caption("No cluster data, re-run the analysis to generate clusters.")

            _phy_resistance = _render_amr_tree_highlight(_phy_active_proj, upload_names, "phy_tree")

            tree_json = newick_to_tree_json(
                newick, clusters=clusters, upload_names=upload_names, resistance=_phy_resistance,
            )
            tree_json_str = json.dumps(tree_json)
            has_clusters = bool(clusters)
            legend_title = "Clusters" if has_clusters else "Sample type"

            html_code = f"""<!DOCTYPE html>
<html>
<head>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  body {{ margin: 0; background: white; font-family: Arial, sans-serif; overflow-x: auto; }}
  .branch {{ fill: none; stroke: #444; stroke-width: 1px; }}
  .leaf-label {{ font-size: 11px; fill: #111; dominant-baseline: middle; }}
  #legend {{
    position: fixed; top: 12px; right: 12px;
    background: rgba(255,255,255,0.96); border: 1px solid #ccc;
    border-radius: 6px; padding: 8px 12px; font-size: 12px; z-index: 10;
  }}
  .legend-row {{ display: flex; align-items: center; gap: 6px; margin: 3px 0; }}
  .legend-dot {{ width: 10px; height: 10px; border-radius: 50%; flex-shrink: 0; }}
</style>
</head>
<body>
<div id="legend"><strong>{legend_title}</strong></div>
<svg id="tree"></svg>
<script>
const data = {tree_json_str};

const labelWidth = 280;
const annotationWidth = 100;
const margin = {{ top: 20, right: labelWidth + annotationWidth, bottom: 30, left: 50 }};
const W = (window.innerWidth || 900);

const root = d3.hierarchy(data);
const nLeaves = root.leaves().length;
const rowH = 20;
const svgH = Math.max(500, nLeaves * rowH) + margin.top + margin.bottom + 40;
const treeW = W - margin.left - margin.right;
const treeH = svgH - margin.top - margin.bottom - 40;

const svg = d3.select('#tree')
    .attr('width', W)
    .attr('height', svgH);

const g = svg.append('g').attr('transform', `translate(${{margin.left}},${{margin.top}})`);

svg.call(d3.zoom()
    .scaleExtent([0.2, 8])
    .on('zoom', e => g.attr('transform', e.transform)));

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

// Rectangular elbow branches
g.selectAll('.branch')
    .data(root.links())
    .join('path')
    .attr('class', 'branch')
    .attr('d', d => `M${{d.source.y}},${{d.source.x}}
                     L${{d.source.y}},${{d.target.x}}
                     L${{d.target.y}},${{d.target.x}}`);

// Internal nodes
g.selectAll('.inode')
    .data(root.descendants().filter(d => d.children))
    .join('circle')
    .attr('cx', d => d.y)
    .attr('cy', d => d.x)
    .attr('r', 2)
    .attr('fill', '#aaa');

// Leaf nodes + labels
const leaves = g.selectAll('.leaf')
    .data(root.leaves())
    .join('g')
    .attr('class', 'leaf')
    .attr('transform', d => `translate(${{d.y}},${{d.x}})`);

const resColor = {{ resistant: '#dc2626', susceptible: '#cbd5e1', no_data: '#e5e7eb' }};
leaves.append('circle')
    .attr('r', d => d.data.resistance === 'resistant' ? 8 : (d.data.resistance ? 4 : 5))
    .attr('fill', d => d.data.resistance ? (resColor[d.data.resistance] || d.data.color) : d.data.color)
    .attr('stroke', d => d.data.resistance === 'resistant' ? '#7f1d1d'
        : (d.data.type === 'leaf_upload' ? '#333' : 'none'))
    .attr('stroke-width', d => d.data.resistance === 'resistant' ? 2 : 1.5)
    .attr('opacity', d => (d.data.resistance && d.data.resistance !== 'resistant') ? 0.45 : 1);

leaves.append('text')
    .attr('class', 'leaf-label')
    .attr('x', 9)
    .text(d => d.data.name);

const clusterLeaves = root.leaves().filter(
    l => l.data.cluster !== undefined &&
         l.data.cluster !== null &&
         String(l.data.cluster) !== ''
);

const clusterIds = Array.from(
    new Set(clusterLeaves.map(l => String(l.data.cluster)))
).sort((a, b) => Number(a) - Number(b));

const clusterColors = {{}};

clusterIds.forEach((c, i) => {{
    clusterColors[c] = d3.hsl(
        (i * 360 / clusterIds.length) % 360,
        0.55,
        0.6
    ).formatHex();
}});

const stripX = treeW + labelWidth;

g.selectAll('.cluster-strip')
    .data(clusterLeaves)
    .join('rect')
    .attr('class', 'cluster-strip')
    .attr('x', stripX)
    .attr('y', d => d.x - 7)
    .attr('width', 14)
    .attr('height', 14)
    .attr('rx', 3)
    .attr('fill', d => clusterColors[String(d.data.cluster)])
    .attr('stroke', '#555')
    .attr('stroke-width', 0.5)
    .append('title')
    .text(d => 'Cluster ' + d.data.cluster);

g.selectAll('.cluster-num')
    .data(clusterLeaves)
    .join('text')
    .attr('class', 'cluster-num')
    .attr('x', stripX + 19)
    .attr('y', d => d.x)
    .attr('dominant-baseline', 'middle')
    .attr('font-size', '10px')
    .attr('fill', '#374151')
    .text(d => d.data.cluster);

if (clusterIds.length) {{
    g.append('text')
        .attr('x', stripX + 7)
        .attr('y', -8)
        .attr('text-anchor', 'middle')
        .attr('font-size', '11px')
        .attr('font-weight', '700')
        .attr('fill', '#374151')
        .text('Cluster');
}}

const legend = document.getElementById('legend');
const typesSeen = {{}};
root.leaves().forEach(n => {{
    const key = n.data.type;
    const label = key === 'leaf_upload' ? 'Upload'
        : key === 'outgroup' ? 'Outgroup'
        : 'Reference';
    if (!typesSeen[key]) {{
        typesSeen[key] = true;
        const row = document.createElement('div');
        row.className = 'legend-row';
        row.innerHTML = `<span class="legend-dot" style="background:${{n.data.color}}"></span>${{label}}`;
        legend.appendChild(row);
    }}
}});
if (clusterIds.length) {{
    const sepC = document.createElement('div');
    sepC.style.cssText = 'border-top:1px solid #ddd;margin:6px 0;';
    legend.appendChild(sepC);
    clusterIds.forEach(c => {{
        const row = document.createElement('div');
        row.className = 'legend-row';
        row.innerHTML = `<span class="legend-dot" style="background:${{clusterColors[c]}};border-radius:3px;"></span>Cluster ${{c}}`;
        legend.appendChild(row);
    }});
}}
const resistanceSeen = new Set();
root.leaves().forEach(n => {{ if (n.data.resistance) resistanceSeen.add(n.data.resistance); }});
if (resistanceSeen.size) {{
    const resLabels = {{ resistant: 'Resistant', susceptible: 'Susceptible', no_data: 'No AMR data' }};
    const sep = document.createElement('div');
    sep.style.cssText = 'border-top:1px solid #ddd;margin:6px 0;';
    legend.appendChild(sep);
    ['resistant', 'susceptible', 'no_data'].forEach(k => {{
        if (!resistanceSeen.has(k)) return;
        const row = document.createElement('div');
        row.className = 'legend-row';
        row.innerHTML = `<span class="legend-dot" style="background:${{resColor[k]}}"></span>${{resLabels[k]}}`;
        legend.appendChild(row);
    }});
}}
</script>
</body>
</html>"""

            st.components.v1.html(html_code, height=950, scrolling=True)

        

            st.download_button(
                label="Download tree (Newick)",
                data=newick,
                file_name="tree.nwk",
                mime="text/plain",
            )

            with st.expander("Raw Newick string"):
                st.text_area("Newick", newick, height=150, label_visibility="collapsed")

    cluster_report = st.session_state.get("phy_cluster_report", {})
    matrix_path    = st.session_state.get("phy_matrix_path")

    if cluster_report:
        st.divider()
        st.subheader("Cluster Report")

        _has_typing  = any("typing" in v for v in cluster_report.values())
        _has_cgmlst  = any("cgmlst" in v for v in cluster_report.values())

        if _has_cgmlst:
            st.caption(
                "Primary clustering: cgMLST allele-based single-linkage at ≤400 AD "
                "(genogroup), a fixed threshold. "
                "SNP nearest-neighbour distance shown as complementary context."
            )
        else:
            st.caption(
                "No cgMLST genogroup available for this run."
            )

        _cr_rows = []

        for _sname, _info in sorted(cluster_report.items()):
            _typing  = _info.get("typing", {})
            _cg      = _info.get("cgmlst", {})
            _cg_geo  = _cg.get("genogroup")
            _row = {
                "Sample": _sname,
                "Type":   "Uploaded" if _info.get("is_uploaded") else "Reference",
                "Input":  _info.get("input_type", "assembly"),
            }
            if _has_cgmlst:
                _row["Cluster"] = str(_cg_geo) if _cg_geo is not None else "—"
            _row["Nearest neighbour"]  = _info.get("nearest_neighbour", "—")
            _row["NN distance (SNPs)"] = f"{_info.get('nn_distance', 0.0):.0f}"
            _row["NN distance (%)"]    = f"{_info.get('nn_distance_pct', 0.0):.4f}%"
            _row["NN type"]            = "Uploaded" if _info.get("nn_is_uploaded") else "Reference"
            if _has_typing:
                _row["MLST ST"]     = _typing.get("mlst_st", "—")
                _row["NG-STAR ST"]  = _typing.get("ngstar_st", "—")
                _row["Mosaic penA"] = _typing.get("mosaic_pena", "—")
                _row["NG-MAST ST"]  = _typing.get("ngmast_st", "—")
            _cr_rows.append(_row)

        _df_cr = pd.DataFrame(_cr_rows)

        def _cr_row_style(row):
            if row["Type"] == "Uploaded":
                return ["background-color:#fef9c3;color:#713f12"] * len(row)
            return [""] * len(row)

        st.dataframe(
            _df_cr.style.apply(_cr_row_style, axis=1),
            use_container_width=True,
            hide_index=True,
        )

        _cr_dl1, _cr_dl2 = st.columns([2, 2])
        _cr_dl1.download_button(
            "⬇ Cluster report (CSV)",
            data=_df_cr.to_csv(index=False).encode(),
            file_name="cluster_report.csv",
            mime="text/csv",
            key="cr_csv_dl",
        )
        if matrix_path and Path(matrix_path).exists():
            _cr_dl2.download_button(
                "⬇ Distance matrix (CSV)",
                data=Path(matrix_path).read_bytes(),
                file_name="distance_matrix.csv",
                mime="text/csv",
                key="cr_matrix_dl",
            )

        _uploaded_rows = [r for r in _cr_rows if r["Type"] == "Uploaded"]
        if _uploaded_rows:
            st.markdown("**Uploaded samples, cluster assignment:**")
            _uc_rows = []
            for _r in _uploaded_rows:
                _uc = {"Sample": _r["Sample"]}
                if _has_cgmlst:
                    _uc["Cluster"] = _r.get("Cluster", "—")
                _uc["Nearest reference"] = _r["Nearest neighbour"] if _r["NN type"] == "Reference" else "—"
                _uc["Distance (SNPs)"]   = _r["NN distance (SNPs)"]
                _uc["Distance (%)"]      = _r["NN distance (%)"]
                if _has_typing:
                    _uc["MLST ST"]     = _r.get("MLST ST", "—")
                    _uc["NG-STAR ST"]  = _r.get("NG-STAR ST", "—")
                    _uc["Mosaic penA"] = _r.get("Mosaic penA", "—")
                    _uc["NG-MAST ST"]  = _r.get("NG-MAST ST", "—")
                _uc_rows.append(_uc)
            st.dataframe(pd.DataFrame(_uc_rows), use_container_width=True, hide_index=True)

    _render_cgmlst_clustering(st.session_state.get("phy_cgmlst_result"), key_prefix="phy_cg")


