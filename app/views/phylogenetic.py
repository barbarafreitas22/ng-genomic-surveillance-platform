from __future__ import annotations
import json
import os
import re
import tempfile
import base64
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

import altair as alt
import pandas as pd
import streamlit as st

from views.shared import *
from views.shared import _inline_metadata_widget


def render() -> None:

    st.html("""
<div style="padding:0.6rem 0 1rem 0;">
  <div style="display:inline-flex;align-items:center;background:rgba(0,0,0,0.04);border:1px solid rgba(0,0,0,0.12);
  border-radius:20px;padding:0.22rem 0.85rem;font-size:0.7rem;font-weight:700;color:#64748b;
  letter-spacing:0.12em;text-transform:uppercase;margin-bottom:0.6rem;">Module 3</div>
  <div style="font-size:1.6rem;font-weight:800;color:#0f172a;">3. Phylogenetic Analysis</div>
</div>
""")

    from backend.phylogeny.viewer import load_tree_newick, newick_to_tree_json


    st.markdown("""
    This module builds whole-genome phylogenetic trees for *Neisseria gonorrhoeae* isolates,
    placing them in context with a curated backbone panel of international reference genomes
    and Portuguese clinical isolates.

    **Fast NJ (Neighbour-Joining)**, uses assembled genomes (FASTA) or raw paired-end FASTQ
    reads. FASTA genomes go straight to sketching; FASTQ reads are optionally trimmed with
    fastp and then sketched directly with SKA2, no assembly step. Faster, uses reads directly,
    for guaranteed accuracy with full QC and assembly, use the Full Pipeline instead. Computes
    pairwise SNP distances with SKA2 and builds a Neighbour-Joining tree with RapidNJ.
    Single-linkage clustering on the SNP distance matrix is applied at three thresholds:
    **≤20 SNPs** (transmission pair), **≤200 SNPs** (outbreak), and **≤2000 SNPs** (stable
    genogroup).

    Sequence typing (MLST · NG-STAR · NG-MAST) and allele-based **cgMLST clustering**
    (chewBBACA) run automatically for assembled-genome samples, applying two epidemiologically
    calibrated thresholds: **7 AD** (outbreak / transmission pair) and **400 AD** (stable
    genogroup). These require an assembly, so samples uploaded as raw reads are shown without
    typing/cgMLST data. chewBBACA and the *N. gonorrhoeae* cgMLST schema are installed and
    downloaded automatically during the Docker image build, no setup is required on the host
    machine.
    """)

    st.caption(
        "Upload assembled genomes (FASTA) and/or paired-end FASTQ reads, in any combination. "
        "SKA2 computes pairwise SNP distances between all genomes in seconds; RapidNJ builds "
        "the Neighbour-Joining tree. Genogroup clusters are assigned by nearest-neighbour SNP "
        "distance against the backbone panel."
    )
    uploaded_genomes = st.file_uploader(
        "Upload assembled genomes (.fa / .fasta / .fna)",
        accept_multiple_files=True,
        type=["fa", "fasta", "fna"],
        key="phy_fasta_upload",
    )
    uploaded_reads = st.file_uploader(
        "Upload paired-end FASTQ reads (.fastq / .fq / .fastq.gz / .fq.gz)",
        accept_multiple_files=True,
        type=["fastq", "fq", "gz"],
        key="phy_fastq_upload",
    )

    _phy_read_pairs: dict = {}
    _phy_trim_reads = True
    if uploaded_reads:
        _phy_read_pairs = group_paired_end(uploaded_reads)
        _phy_pe_only = {sid: g for sid, g in _phy_read_pairs.items() if g.get("R1") and g.get("R2")}
        _n_unpaired = len(_phy_read_pairs) - len(_phy_pe_only)
        if _n_unpaired:
            st.warning(f"{_n_unpaired} sample(s) could not be paired (R1/R2 not detected) — excluded.")
        _phy_read_pairs = _phy_pe_only
        st.caption(f"{len(uploaded_reads)} file(s) → {len(_phy_read_pairs)} paired-end sample(s) detected.")
        _phy_trim_reads = st.checkbox(
            "Trim reads with fastp before sketching (recommended)",
            value=True,
            key="phy_trim_reads",
            help="Removes adapters and low-quality bases before SKA2 sketching. Untrimmed "
                 "adapter sequence can introduce false k-mers into the distance matrix. "
                 "Uncheck to sketch raw reads directly (faster, less safe).",
        )

    _phy_all_sids = [Path(f.name).stem for f in uploaded_genomes] + list(_phy_read_pairs.keys())
    if _phy_all_sids:
        _phy_proj = st.session_state.get("active_project")
        _inline_metadata_widget(_phy_all_sids, _phy_proj, "phy")

    _phy_has_input = bool(uploaded_genomes) or bool(_phy_read_pairs)
    _phy_busy = bool(st.session_state.get("_phy_pending_job"))
    if _phy_has_input and st.button("Run Fast NJ", disabled=_phy_busy):
        from backend.paths import RUNS_DIR as _RUNS_DIR
        import uuid as _uuid
        _stage = Path(_RUNS_DIR) / f"stage_{_uuid.uuid4().hex}"
        _stage.mkdir(parents=True, exist_ok=True)

        _samples_payload = []
        for _f in uploaded_genomes:
            _dest = _stage / _f.name
            _dest.write_bytes(_f.read())
            _samples_payload.append({"kind": "fasta", "path": str(_dest)})

        for _sid, _pair in _phy_read_pairs.items():
            _r1p = _stage / _pair["R1"].name
            _r2p = _stage / _pair["R2"].name
            _r1p.write_bytes(_pair["R1"].getbuffer())
            _r2p.write_bytes(_pair["R2"].getbuffer())
            _samples_payload.append({
                "kind": "fastq", "name": sanitize_name(_sid),
                "r1": str(_r1p), "r2": str(_r2p), "trim": _phy_trim_reads,
            })

        _jid = db.submit_job(
            "_phylogeny", "_phy_nj", "phylogeny_nj", {"samples": _samples_payload}
        )
        st.session_state["_phy_pending_job"] = _jid
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
                st.info("Queued, waiting for worker…")
            else:
                with st.status("Running phylogenetic analysis…", state="running", expanded=True):
                    st.write("Building the SKA2 tree and computing cgMLST allele-based clustering.")
                    st.write(
                        "The cgMLST step (chewBBACA allele calling) is the slowest part and can take "
                        "several minutes for larger batches — this page updates automatically, no need "
                        "to resubmit."
                    )
        elif _phy_jstatus == "done" and _phy_job.get("result"):
            _phy_res  = json.loads(_phy_job["result"])
            st.session_state["phy_run_id"]         = _phy_res["run_id"]
            st.session_state["phy_tree_path"]      = _phy_res["tree_path"]
            st.session_state["phy_tree_method"]    = "Fast NJ"
            st.session_state["phy_clusters"]       = _phy_res.get("clusters", {})
            st.session_state["phy_cluster_report"] = _phy_res.get("cluster_report", {})
            st.session_state["phy_matrix_path"]    = _phy_res.get("matrix_path")
            st.session_state["phy_upload_names"]   = _phy_res.get("upload_names", [])
            st.session_state["phy_cgmlst_result"]  = _phy_res.get("cgmlst_result")
            st.session_state["phy_snp_clusters"]   = _phy_res.get("snp_clusters")
            st.session_state.pop("_phy_pending_job", None)
            st.rerun()
        elif _phy_jstatus == "error":
            st.error(f"Phylogeny error: {_phy_job.get('error_msg', 'Unknown error')}")
            st.session_state.pop("_phy_pending_job", None)

    if st.session_state.get("_phy_pending_job") and db:
        _phy_poll_status()

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
            _gg_clusters = _cgmlst_for_tree.get("genogroup_clusters") or {}
            if _gg_clusters:
                clusters = _gg_clusters
                cluster_source = "cgMLST genogroup (≤400 AD)"
            else:
                clusters = st.session_state.get("phy_clusters", {})
                cluster_source = "SNP (ska2)"
            upload_names = st.session_state.get("phy_upload_names", [])

            if clusters:
                from collections import Counter
                ctr = Counter(clusters.values())
                n_multi = sum(1 for v in ctr.values() if v >= 2)
                st.caption(f"{len(ctr)} clusters detected ({cluster_source}), {n_multi} with ≥2 members")
            else:
                st.caption("No cluster data, re-run the analysis to generate clusters.")

            tree_json = newick_to_tree_json(newick, clusters=clusters, upload_names=upload_names)
            tree_json_str = json.dumps(tree_json)
            has_clusters = bool(clusters)
            legend_title = "Clusters" if has_clusters else "Sample type"

            html_code = f"""<!DOCTYPE html>
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
</style>
</head>
<body>
<div id="legend"><strong>{legend_title}</strong></div>
<svg id="tree"></svg>
<script>
const data = {tree_json_str};

const margin = {{ top: 20, right: 460, bottom: 30, left: 50 }};
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

// Cluster shadows
const clusterMap = {{}};
root.leaves().forEach(l => {{
    const c = l.data.cluster;
    if (c && c !== '') {{
        if (!clusterMap[c]) clusterMap[c] = [];
        clusterMap[c].push(l);
    }}
}});

const shadowPalette = [
    'rgba(147,197,253,0.30)', 'rgba(249,168,212,0.30)',
    'rgba(134,239,172,0.30)', 'rgba(253,224,132,0.30)',
    'rgba(196,181,253,0.30)', 'rgba(103,232,249,0.30)',
    'rgba(254,202,202,0.30)', 'rgba(167,243,208,0.30)',
];
let ci = 0;
Object.entries(clusterMap)
    .filter(([, ls]) => ls.length >= 2)
    .sort(([a],[b]) => +a - +b)
    .forEach(([c, ls]) => {{
        const xs = ls.map(l => l.x);
        const yMin = Math.min(...xs) - rowH * 0.6;
        const yMax = Math.max(...xs) + rowH * 0.6;
        const col = shadowPalette[ci++ % shadowPalette.length];
        g.append('rect')
            .attr('x', 0)
            .attr('y', yMin)
            .attr('width', treeW + 10)
            .attr('height', yMax - yMin)
            .attr('rx', 8)
            .attr('fill', col);
        const labelX = treeW + 240;
        const labelY = (yMin + yMax) / 2;
        const labelH = 22;
        g.append('rect')
            .attr('x', labelX - 8)
            .attr('y', labelY - labelH / 2)
            .attr('width', 110)
            .attr('height', labelH)
            .attr('rx', 5)
            .attr('fill', col.replace('0.30', '0.70'))
            .attr('stroke', col.replace('0.30', '0.90'))
            .attr('stroke-width', 1);
        g.append('text')
            .attr('x', labelX)
            .attr('y', labelY)
            .attr('dominant-baseline', 'middle')
            .attr('font-size', '12px')
            .attr('font-weight', '700')
            .attr('fill', '#1f2937')
            .text('Cluster ' + c);
    }});

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

leaves.append('circle')
    .attr('r', 5)
    .attr('fill', d => d.data.color)
    .attr('stroke', d => d.data.type === 'leaf_upload' ? '#333' : 'none')
    .attr('stroke-width', 1.5);

leaves.append('text')
    .attr('class', 'leaf-label')
    .attr('x', 9)
    .text(d => d.data.name);

// Legend — node types (clusters are labelled on the shadows)
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
</script>
</body>
</html>"""

            st.components.v1.html(html_code, height=950, scrolling=True)

        

            st.download_button(
                label="⬇ Download tree (Newick)",
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
                "Primary clustering: cgMLST allele-based single-linkage at ≤7 AD (outbreak) "
                "and ≤400 AD (genogroup), fixed, dataset-independent thresholds. "
                "SNP nearest-neighbour distance shown as complementary context."
            )
        else:
            st.caption(
                "Clustering: single-linkage on SKA2 pairwise SNP distances. "
                "Nearest-neighbour distance shown for each sample against the full panel."
            )

        _cr_rows = []
        _cg_outbreak_counts: dict = {}
        if _has_cgmlst:
            for _v in cluster_report.values():
                _oc = (_v.get("cgmlst") or {}).get("outbreak_cluster")
                if _oc is not None:
                    _cg_outbreak_counts[_oc] = _cg_outbreak_counts.get(_oc, 0) + 1
        _outbreak_sids: set = set()

        for _sname, _info in sorted(cluster_report.items()):
            _typing  = _info.get("typing", {})
            _cg      = _info.get("cgmlst", {})
            _cg_out  = _cg.get("outbreak_cluster")
            _cg_geo  = _cg.get("genogroup")
            if _cg_out is not None and _cg_outbreak_counts.get(_cg_out, 0) >= 2:
                _outbreak_sids.add(_sname)
            _row = {
                "Sample": _sname,
                "Type":   "Uploaded" if _info.get("is_uploaded") else "Reference",
                "Input":  _info.get("input_type", "assembly"),
            }
            if _has_cgmlst:
                _row["cgMLST Genogroup (400 AD)"] = str(_cg_geo) if _cg_geo is not None else "—"
                _row["cgMLST Outbreak (7 AD)"]    = str(_cg_out) if _cg_out is not None else "—"
            _row["Nearest neighbour"]  = _info.get("nearest_neighbour", "—")
            _row["NN distance (SNPs)"] = f"{_info.get('nn_distance', 0.0):.0f}"
            _row["NN distance (%)"]    = f"{_info.get('nn_distance_pct', 0.0):.4f}%"
            _row["NN type"]            = "Uploaded" if _info.get("nn_is_uploaded") else "Reference"
            if not _has_cgmlst:
                _row["SNP Genogroup"] = str(_info.get("cluster", "—"))
            if _has_typing:
                _row["MLST ST"]           = _typing.get("mlst_st", "—")
                _row["NG-STAR CC"]        = _typing.get("ngstar_cc", "—")
                _row["NG-MAST Genogroup"] = _typing.get("ngmast_genogroup", "—")
                _row["Mosaic penA"]       = _typing.get("mosaic_pena", "—")
            _cr_rows.append(_row)

        _df_cr = pd.DataFrame(_cr_rows)

        def _cr_row_style(row):
            if row["Sample"] in _outbreak_sids:
                return ["background-color:#fef2f2;color:#7f1d1d"] * len(row)
            if row["Type"] == "Uploaded":
                return ["background-color:#fef9c3;color:#713f12"] * len(row)
            return [""] * len(row)

        if _outbreak_sids:
            st.error(
                f"{len(_outbreak_sids)} uploaded sample(s) share ≤7 allele differences "
                "(cgMLST outbreak signal), verify with epidemiological data."
            )

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
                    _uc["cgMLST Genogroup (400 AD)"] = _r.get("cgMLST Genogroup (400 AD)", "—")
                    _uc["cgMLST Outbreak (7 AD)"]    = _r.get("cgMLST Outbreak (7 AD)", "—")
                else:
                    _uc["SNP Genogroup"] = _r.get("SNP Genogroup", "—")
                _uc["Nearest reference"] = _r["Nearest neighbour"] if _r["NN type"] == "Reference" else "—"
                _uc["Distance (SNPs)"]   = _r["NN distance (SNPs)"]
                _uc["Distance (%)"]      = _r["NN distance (%)"]
                if _has_typing:
                    _uc["MLST ST"]           = _r.get("MLST ST", "—")
                    _uc["NG-STAR CC"]        = _r.get("NG-STAR CC", "—")
                    _uc["NG-MAST Genogroup"] = _r.get("NG-MAST Genogroup", "—")
                    _uc["Mosaic penA"]       = _r.get("Mosaic penA", "—")
                _uc_rows.append(_uc)
            st.dataframe(pd.DataFrame(_uc_rows), use_container_width=True, hide_index=True)

    _cgmlst_result = st.session_state.get("phy_cgmlst_result")
    if _cgmlst_result:
        st.divider()
        st.subheader("cgMLST Clustering")
        st.caption(
            f"Core genome: {_cgmlst_result.get('core_loci_count', '?')} loci (≥95% presence). "
            "Single-linkage clustering applied at two calibrated thresholds. "
            "Thresholds are fixed and dataset-independent, enabling longitudinal comparison."
        )

        _cg_stats = _cgmlst_result.get("sample_stats", {})
        _cg_outbreak = _cgmlst_result.get("outbreak_clusters", {})
        _cg_genogroup = _cgmlst_result.get("genogroup_clusters", {})

        if _cg_stats:
            _cg_tab1, _cg_tab2 = st.tabs(
                ["Outbreak Clusters (≤7 AD)", "Genogroup Assignment (≤400 AD)"]
            )

            with _cg_tab1:
                st.caption(
                    "**7 allele difference threshold**, captures likely transmission pairs "
                    "and active outbreaks. Samples sharing ≤7 AD are co-clustered. "
                    "Two samples are in the same outbreak cluster if their cgMLST distance "
                    "is ≤7 AD, or transitively linked through intermediates."
                )
                _oc_rows = []
                for _sid, _st in sorted(_cg_stats.items()):
                    _oc_rows.append({
                        "Sample":           _sid,
                        "Outbreak cluster": _st.get("outbreak_cluster", "—"),
                        "Loci assigned":    _st.get("loci_assigned", "—"),
                        "Loci total":       _st.get("loci_total", "—"),
                        "% assigned":       f"{_st.get('pct_assigned', 0):.1f}%",
                    })
                _df_oc = pd.DataFrame(_oc_rows)

                _oc_cluster_sizes = {}
                for _sid, _cl in _cg_outbreak.items():
                    _oc_cluster_sizes[_cl] = _oc_cluster_sizes.get(_cl, 0) + 1
                _outbreak_sids = {s for s, c in _cg_outbreak.items() if _oc_cluster_sizes.get(c, 1) > 1}

                def _oc_style(row):
                    return (
                        ["background-color:#fef2f2;color:#7f1d1d"] * len(row)
                        if row["Sample"] in _outbreak_sids
                        else [""] * len(row)
                    )

                st.dataframe(
                    _df_oc.style.apply(_oc_style, axis=1),
                    use_container_width=True,
                    hide_index=True,
                )

                _outbreak_clusters = {c for c, n in _oc_cluster_sizes.items() if n > 1}
                if _outbreak_clusters:
                    st.warning(
                        f"{len(_outbreak_sids)} sample(s) share ≤7 allele differences, "
                        "possible active transmission. Verify with epidemiological data."
                    )
                else:
                    st.success("No samples within 7 allele differences, no outbreak signal detected.")

            with _cg_tab2:
                st.caption(
                    "**400 allele difference threshold (Ng_cgc_400)**, stable genogroup "
                    "assignment persisting over decades. Equivalent to lineage-level classification. "
                    "Use for longitudinal surveillance and strain tracking."
                )
                _gg_rows = []
                _gg_cluster_sizes = {}
                for _sid, _cl in _cg_genogroup.items():
                    _gg_cluster_sizes[_cl] = _gg_cluster_sizes.get(_cl, 0) + 1
                for _sid, _st in _cg_stats.items():
                    _gg_rows.append({
                        "Sample":    _sid,
                        "Genogroup": _st.get("genogroup"),
                        "% assigned": f"{_st.get('pct_assigned', 0):.1f}%",
                    })
                _gg_rows.sort(key=lambda r: (
                    r["Genogroup"] if r["Genogroup"] is not None else float("inf"),
                    r["Sample"],
                ))
                for _row in _gg_rows:
                    _row["Genogroup"] = _row["Genogroup"] if _row["Genogroup"] is not None else "—"
                st.dataframe(pd.DataFrame(_gg_rows), use_container_width=True, hide_index=True)

            _cg_pairs = _cgmlst_result.get("pairwise", [])
            if _cg_pairs:
                with st.expander("Pairwise allele differences"):
                    _pair_rows = []
                    for _p in sorted(_cg_pairs, key=lambda x: x["allele_diff"]):
                        _pair_rows.append({
                            "Sample A":          _p["sample_a"],
                            "Sample B":          _p["sample_b"],
                            "Allele differences": _p["allele_diff"],
                            "Outbreak link (≤7)": "Yes" if _p["outbreak_link"] else "No",
                            "Same genogroup (≤400)": "Yes" if _p["same_genogroup"] else "No",
                        })
                    _df_pairs = pd.DataFrame(_pair_rows)

                    def _pair_style(row):
                        if row["Outbreak link (≤7)"] == "Yes":
                            return ["background-color:#fef2f2;color:#7f1d1d"] * len(row)
                        return [""] * len(row)

                    st.dataframe(
                        _df_pairs.style.apply(_pair_style, axis=1),
                        use_container_width=True,
                        hide_index=True,
                    )

            _cg_matrix_csv = _cgmlst_result.get("matrix_csv")
            if _cg_matrix_csv and Path(_cg_matrix_csv).exists():
                st.download_button(
                    "⬇ cgMLST distance matrix (CSV)",
                    data=Path(_cg_matrix_csv).read_bytes(),
                    file_name="cgmlst_distances.csv",
                    mime="text/csv",
                    key="cgmlst_matrix_dl",
                )

    _snp_clusters = st.session_state.get("phy_snp_clusters")
    if _snp_clusters:
        st.markdown("---")
        st.markdown("#### WGS SNP Clustering (ska2)")
        st.caption(
            f"Single-linkage clustering on pairwise split k-mer SNP distances. "
            f"Transmission threshold: ≤{_snp_clusters.get('transmission', {}).get('threshold', SNP_TRANSMISSION_THRESHOLD)} SNPs · "
            f"Outbreak threshold: ≤{_snp_clusters.get('outbreak', {}).get('threshold', SNP_OUTBREAK_THRESHOLD)} SNPs."
        )

        _snp_tabs = st.tabs([
            f"Transmission clusters (≤{_snp_clusters.get('transmission', {}).get('threshold', SNP_TRANSMISSION_THRESHOLD)} SNPs)",
            f"Outbreak clusters (≤{_snp_clusters.get('outbreak', {}).get('threshold', SNP_OUTBREAK_THRESHOLD)} SNPs)",
        ])

        for _snp_tab, _snp_level in zip(_snp_tabs, ["transmission", "outbreak"]):
            with _snp_tab:
                _level_data = _snp_clusters.get(_snp_level, {})
                _level_clusters = _level_data.get("clusters", {})
                _level_pairs    = _level_data.get("pairs", {})

                if not _level_clusters:
                    st.info("No data available.")
                else:
                    from collections import defaultdict as _ddict
                    _by_cluster: dict = _ddict(list)
                    for _samp, _clust in sorted(_level_clusters.items()):
                        _by_cluster[_clust].append(_samp)

                    _multi = {k: v for k, v in _by_cluster.items() if len(v) > 1}
                    _single = {k: v for k, v in _by_cluster.items() if len(v) == 1}

                    if _multi:
                        st.markdown(f"**{len(_multi)} cluster(s) with ≥2 samples:**")
                        for _cid, _members in sorted(_multi.items()):
                            with st.expander(f"Cluster {_cid}, {len(_members)} samples", expanded=True):
                                _pair_rows = []
                                for _i, _a in enumerate(_members):
                                    for _b in _members[_i + 1:]:
                                        _key = f"{_a}||{_b}"
                                        _rkey = f"{_b}||{_a}"
                                        _dist = _level_pairs.get(_key) or _level_pairs.get(_rkey)
                                        if _dist is not None:
                                            _pair_rows.append({"Sample A": _a, "Sample B": _b, "SNPs": _dist})
                                if _pair_rows:
                                    st.dataframe(
                                        pd.DataFrame(_pair_rows).sort_values("SNPs"),
                                        hide_index=True,
                                        use_container_width=True,
                                    )
                                else:
                                    st.write(", ".join(_members))
                    else:
                        st.info("No samples clustered together at this threshold, all are genetically distinct.")

                    if _single:
                        st.markdown(f"**{len(_single)} singleton(s):** " + ", ".join(
                            v[0] for v in sorted(_single.values())
                        ))


