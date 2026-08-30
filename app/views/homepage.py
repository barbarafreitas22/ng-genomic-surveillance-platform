from __future__ import annotations
import json
import os
import time
from pathlib import Path

import pandas as pd
import streamlit as st

from views.shared import *
from views.shared import (
    _load_project_cached, _rebuild_pipeline_results,
    _inline_metadata_widget, _clear_project_results,
)
from backend.models import AMRResult as _AMRResult


def render() -> None:
    st.markdown("""
<div style="
  background: #f8fafc;
  border: 1px solid #e2e8f0;
  border-radius: 16px;
  padding: 2.4rem 2.8rem 2rem 2.8rem;
  margin-bottom: 1.4rem;
">
  <div style="display:flex;gap:0.5rem;flex-wrap:wrap;margin-bottom:0.9rem;">
    <div style="display:inline-block;background:rgba(0,0,0,0.04);border:1px solid rgba(0,0,0,0.12);
      border-radius:20px;padding:0.22rem 0.85rem;font-size:0.7rem;font-weight:700;color:#64748b;
      letter-spacing:0.12em;text-transform:uppercase;">
      Whole Genome Sequencing
    </div>
    <div style="display:inline-block;background:rgba(0,0,0,0.04);border:1px solid rgba(0,0,0,0.14);
      border-radius:20px;padding:0.22rem 0.85rem;font-size:0.7rem;font-weight:700;color:#64748b;
      letter-spacing:0.12em;text-transform:uppercase;">
      AMR Surveillance
    </div>
  </div>

  <div style="font-size:2.4rem;font-weight:800;color:#0f172a;margin:0 0 0.25rem 0;
    line-height:1.15;letter-spacing:-0.03em;">
    NGsurv
  </div>
  <div style="font-size:1.05rem;font-weight:600;color:#64748b;margin:0 0 0.75rem 0;line-height:1.4;">
    <em style="font-style:italic;">Neisseria gonorrhoeae</em> Genomic Surveillance Platform
  </div>

  <p style="color:#64748b;font-size:0.88rem;margin:0;max-width:780px;line-height:1.7;">
    This platform offers a end-to-end Whole Genome Sequencing analysis. When raw reads are provided,
    quality control and genome assembly are performed before antimicrobial resistance profiling,
    phylogenetic contextualization, and genomic clinical interpretation. Individual modules can be accessed
    via the sidebar, and the complete automated workflow is available below.
  </p>
</div>
""", unsafe_allow_html=True)

    st.html("""
<style>
  .wf-wrap {
    margin: 0 0 1rem 0;
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
  }
  .wf-label {
    font-size: 0.6rem;
    font-weight: 700;
    letter-spacing: 0.15em;
    text-transform: uppercase;
    color: #94a3b8;
    margin-bottom: 0.8rem;
  }
  .wf-row {
    display: flex;
    align-items: stretch;
    gap: 0;
    overflow-x: auto;
    padding-bottom: 0.3rem;
  }
  .wf-step {
    flex-shrink: 0;
    background: rgba(0,0,0,0.03);
    border: 1px solid rgba(0,0,0,0.12);
    border-radius: 12px;
    padding: 0.8rem 1.1rem;
    text-align: center;
    min-width: 108px;
  }
  .wf-step.last {
    border-color: rgba(0,0,0,0.18);
  }
  .wf-title {
    font-size: 0.78rem;
    font-weight: 700;
    color: #334155;
    white-space: nowrap;
  }
  .wf-step.last .wf-title { color: #475569; }
  .wf-sub {
    font-size: 0.62rem;
    color: #94a3b8;
    margin-top: 0.15rem;
    white-space: nowrap;
  }
  .wf-step.first .wf-sub { color: #94a3b8; }
  .wf-arrow {
    flex-shrink: 0;
    display: flex;
    align-items: center;
    padding: 0 0.5rem;
    color: #cbd5e1;
    font-size: 1.2rem;
  }
</style>
<div class="wf-wrap">
  <div class="wf-label">End-to-end automated workflow</div>
  <div class="wf-row">
    <div class="wf-step first">
      <div class="wf-title">Raw FASTQ</div>
      <div class="wf-sub">Sequencer output</div>
    </div>
    <div class="wf-arrow">&#8594;</div>
    <div class="wf-step">
      <div class="wf-title">QC &amp; Cleaning</div>
      <div class="wf-sub">fastp · Kraken2</div>
    </div>
    <div class="wf-arrow">&#8594;</div>
    <div class="wf-step">
      <div class="wf-title">Assembly</div>
      <div class="wf-sub">SPAdes · Biopython</div>
    </div>
    <div class="wf-arrow">&#8594;</div>
    <div class="wf-step">
      <div class="wf-title">AMR Profiling</div>
      <div class="wf-sub">Minimap2 · European 2020</div>
    </div>
    <div class="wf-arrow">&#8594;</div>
    <div class="wf-step">
      <div class="wf-title">Phylogenetics</div>
      <div class="wf-sub">SKA2 · MLST · NG-STAR</div>
    </div>
    <div class="wf-arrow">&#8594;</div>
    <div class="wf-step last">
      <div class="wf-title">Clinical Report</div>
      <div class="wf-sub">AMR · Alerts</div>
    </div>
  </div>
</div>
""")

    st.divider()

    left_col, right_col = st.columns([1, 1], gap="large")

    with left_col:
        st.markdown("""
<div style="font-size:0.7rem;font-weight:700;text-transform:uppercase;letter-spacing:0.1em;
            color:#64748b;margin-bottom:0.5rem;">Project Setup</div>
""", unsafe_allow_html=True)

        project_dirs = list_subdirs(PROJECTS_DIR)
        project_list = [p.name for p in project_dirs]

        mode = st.radio(
            "Mode",
            ["Use existing project", "Create new project"],
            horizontal=True,
            label_visibility="collapsed",
        )

        if mode == "Use existing project":
            selected_project = st.selectbox(
                "Select project",
                project_list if project_list else ["(No projects available)"],
            )
            if project_list:
                _db_data = _load_project_cached(selected_project)
                if selected_project != st.session_state.get("active_project"):
                    _clear_project_results()
                    st.session_state["active_project"] = selected_project
                    if _db_data.get("samples"):
                        st.session_state["full_pipeline_results"] = _rebuild_pipeline_results(_db_data)
                        if _db_data.get("phylogeny"):
                            _phy_db = _db_data["phylogeny"]
                            st.session_state["last_tree_path"]      = _phy_db.get("tree_path")
                            st.session_state["last_cluster_report"] = _phy_db.get("cluster_report", {})
                            st.session_state["last_matrix_path"]    = _phy_db.get("matrix_path")
                            st.session_state["last_cgmlst_result"]  = _phy_db.get("cgmlst_result")
                            _cr_db = _phy_db.get("cluster_report", {})
                            st.session_state["last_clusters"] = {
                                n: i["cluster"] for n, i in _cr_db.items() if "cluster" in i
                            }
                _n_saved = len(_db_data["samples"])
                if _n_saved:
                    st.success(f"Active project: **{selected_project}** · {_n_saved} {plural(_n_saved, 'sample')} saved")
                else:
                    st.success(f"Active project: **{selected_project}**")
                if st.button("Delete all results", key="del_all_home_btn"):
                    if db is not None:
                        db.delete_qc(selected_project)
                        db.delete_assembly(selected_project)
                        db.delete_amr(selected_project)
                        db.delete_phylogeny(selected_project)
                        db.delete_alerts(selected_project)
                        db.delete_samples(selected_project)
                    _clear_project_results()
                    _load_project_cached.clear()
                    st.rerun()
            else:
                st.warning("No projects found. Create one first.")
            active_project = selected_project if project_list else None

        else:
            new_project_name = st.text_input(
                "Project name",
                placeholder="e.g. NG_Surveillance_2026",
            )
            if st.button("Create Project", use_container_width=True):
                cleaned = sanitize_name(new_project_name)
                if cleaned and cleaned != "sample":
                    proj_dir = PROJECTS_DIR / cleaned
                    for sub in ["raw", "results/qc", "results/assembly",
                                "results/amr", "results/phylogeny"]:
                        ensure_dir(proj_dir / sub)
                    if db is not None:
                        db.init_project(cleaned)
                    _clear_project_results()
                    st.session_state["active_project"] = cleaned
                    st.success(f"Project **{cleaned}** created.")
                    st.rerun()
                else:
                    st.error("Enter a valid project name.")
            active_project = new_project_name if new_project_name else None

    with right_col:
        st.markdown("""
<div style="font-size:0.7rem;font-weight:700;text-transform:uppercase;letter-spacing:0.1em;
            color:#64748b;margin-bottom:0.6rem;">Run Full Pipeline</div>
<div style="display:grid;grid-template-columns:1fr 1fr;gap:0.5rem;margin-bottom:0.75rem;">
  <div style="background:rgba(37,99,235,0.08);border:1px solid rgba(37,99,235,0.28);
              border-radius:10px;padding:0.65rem 0.85rem;">
    <div style="font-weight:700;font-size:0.82rem;color:#60a5fa;margin-bottom:0.25rem;">
      Raw reads (FASTQ)
    </div>
    <div style="font-size:0.72rem;color:#94a3b8;">
      <code style="font-size:0.68rem;">.fastq &nbsp;.fq &nbsp;.fastq.gz</code>
    </div>
  </div>
  <div style="background:rgba(124,58,237,0.08);border:1px solid rgba(124,58,237,0.28);
              border-radius:10px;padding:0.65rem 0.85rem;">
    <div style="font-weight:700;font-size:0.82rem;color:#a78bfa;margin-bottom:0.25rem;">
      Assembled genomes (FASTA)
    </div>
    <div style="font-size:0.72rem;color:#94a3b8;">
      <code style="font-size:0.68rem;">.fasta &nbsp;.fa &nbsp;.fna</code>
    </div>
  </div>
</div>
""", unsafe_allow_html=True)

        _tab_upload, _tab_local = st.tabs(["Upload files", "Local folder"])

        uploaded_home: list = []
        _local_path_files: list[Path] = []

        with _tab_upload:
            uploaded_home = st.file_uploader(
                "Upload FASTQ reads or assembled FASTA genomes",
                type=["fastq", "fq", "gz", "fasta", "fa", "fna", "skf"],
                accept_multiple_files=True,
                key="homepage_upload",
            ) or []

        with _tab_local:
            _local_folder = st.text_input(
                "Folder path",
                placeholder="/Volumes/disco/samples",
                key="fp_local_folder",
                help="Paste the full path to a folder containing FASTQ/FASTA files. "
                     "Files are linked directly — no copying needed.",
            )
            if _local_folder:
                _lp = Path(_local_folder.strip())
                if not _lp.is_dir():
                    st.error(f"Folder not found: {_lp}")
                else:
                    _FASTQ_EXTS = {".fastq", ".fq", ".fastq.gz", ".fq.gz",
                                   ".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz"}
                    _local_path_files = [
                        f for f in sorted(_lp.iterdir())
                        if f.is_file() and "".join(f.suffixes).lower() in _FASTQ_EXTS
                    ]
                    if not _local_path_files:
                        st.warning("No FASTQ/FASTA files found in that folder.")

        _active_files: list = _local_path_files if _local_path_files else uploaded_home

        if _active_files:
            _fp_groups = group_paired_end(_active_files)
            _n_files_rdy = len(_active_files)
            _n_fp_groups = len(_fp_groups)
            st.success(
                f"{_n_files_rdy} {plural(_n_files_rdy, 'file')} ready, "
                f"**{_n_fp_groups} {plural(_n_fp_groups, 'sample')}** detected."
            )
            _fp_id_map: dict = {}
            for _fpi, _fpbase in enumerate(_fp_groups):
                _fp_lbl = "  +  ".join(
                    f"`{_fp_groups[_fpbase][r].name}`"
                    for r in ("R1", "R2") if _fp_groups[_fpbase].get(r)
                )
                _fpc1, _fpc2 = st.columns([3, 2])
                _fpc1.markdown(_fp_lbl)
                _fp_inp = _fpc2.text_input(
                    "ID", value=_fpbase, placeholder="e.g. NG-001",
                    key=f"fp_sid_{_fpi}", label_visibility="collapsed",
                )
                _fp_id_map[_fpbase] = sanitize_name(_fp_inp) if _fp_inp.strip() else _fpbase

            _inline_metadata_widget(
                list(_fp_id_map.values()), active_project, "fp"
            )

        _fp_job_pending = bool(st.session_state.get("_fp_pending_job"))

        if st.button("▶  Run Full Pipeline", type="primary", use_container_width=True,
                     disabled=_fp_job_pending):
            if run_full_pipeline is None:
                st.error("Pipeline backend could not be imported.")
                st.stop()
            if not active_project:
                st.error("Select or create a project first.")
                st.stop()
            if not _active_files:
                st.error("Upload files or enter a local folder path.")
                st.stop()
            if not db:
                st.error("Job queue database is not available.")
                st.stop()

            project_input_dir = project_raw_dir(active_project)
            ensure_dir(project_input_dir)
            _fp_groups_run = group_paired_end(_active_files)
            if not _fp_groups_run:
                st.error("No samples detected. Check if the filenames contain _1/_2 or _R1/_R2.")
                st.stop()

            _fp_id_map_run = {
                _fpbase: sanitize_name(
                    st.session_state.get(f"fp_sid_{_fpi}", _fpbase)
                ) or _fpbase
                for _fpi, _fpbase in enumerate(_fp_groups_run)
            }

            _using_local = bool(_local_path_files)
            input_paths: list[Path] = []
            try:
                for _fpbase, _fppair in _fp_groups_run.items():
                    _sid_run = _fp_id_map_run[_fpbase]
                    for _rk in ("R1", "R2"):
                        _ffile = _fppair.get(_rk)
                        if _ffile:
                            _ext  = "".join(Path(_ffile.name).suffixes)
                            _dest = project_input_dir / f"{_sid_run}_{_rk}{_ext}"
                            if _using_local:
                                _src = Path(_ffile) if isinstance(_ffile, Path) else Path(_ffile.name)
                                if _dest.exists() or _dest.is_symlink():
                                    _dest.unlink()
                                _dest.symlink_to(_src.resolve())
                            else:
                                safe_write_bytes(_dest, _ffile.getbuffer())
                            input_paths.append(_dest)
            except Exception as _write_err:
                st.error(f"Failed to prepare input files: {_write_err}")
                st.stop()

            if not input_paths:
                st.error("No files were saved.")
                st.stop()

            _all_sids_run = list(_fp_id_map_run.values())

            _jid = db.submit_job(
                active_project, "_full_pipeline", "full_pipeline",
                {
                    "input_paths": [str(p) for p in input_paths],
                },
            )
            st.session_state["_fp_pending_job"]     = _jid
            st.session_state["_fp_pending_sids"]    = _all_sids_run
            st.session_state["_fp_pending_started"] = time.time()
            st.rerun()

        _fp_pending_job = st.session_state.get("_fp_pending_job")
        if _fp_pending_job and db:
            _fp_jobs_all = {j["id"]: j for j in db.get_project_jobs(active_project)}
            _fp_job      = _fp_jobs_all.get(_fp_pending_job)
            _fp_sids     = st.session_state.get("_fp_pending_sids", [])
            _fp_n_total  = len(_fp_sids)

            if _fp_job is None:
                st.error("Pipeline job not found — it may have been lost when the worker restarted.")
                for _k in ("_fp_pending_job", "_fp_pending_sids", "_fp_pending_started"):
                    st.session_state.pop(_k, None)

            elif _fp_job["status"] in ("queued", "running"):
                _fp_live = db.load_project(active_project).get("samples", {})
                _fp_started = _fp_job.get("started_at") or _fp_job.get("created_at") or ""
                def _fp_fresh(_sid, _key):
                    return (_fp_live.get(_sid, {}).get(_key) or "") >= _fp_started
                _fp_n_done = sum(1 for _sid in _fp_sids if _fp_fresh(_sid, "amr_ran_at"))

                with st.status(
                    f"Running Full Pipeline: {active_project} · {_fp_n_total} {plural(_fp_n_total, 'sample')}",
                    expanded=True,
                ) as _fp_status:
                    _fp_label = "Queued…" if _fp_job["status"] == "queued" else f"{_fp_n_done}/{_fp_n_total} {plural(_fp_n_total, 'sample')}"
                    st.caption(_fp_label)
                    _fp_status.update(label=_fp_label)
                time.sleep(3)
                st.rerun()

            elif _fp_job["status"] == "done" and _fp_job.get("result"):
                results = json.loads(_fp_job["result"])
                results["amr"] = {
                    sid: (_AMRResult.from_dict(v, sample_id=sid)
                          if isinstance(v, dict) and "error" not in v else v)
                    for sid, v in results.get("amr", {}).items()
                }

                _elapsed = time.time() - st.session_state.get("_fp_pending_started", time.time())
                st.caption(f"⏱ Pipeline wall time: **{_elapsed:.1f}s**")

                if results["phylogeny"].get("error"):
                    st.write(f"⚠️ Phylogeny: {results['phylogeny']['error']}")

                st.session_state["full_pipeline_results"] = results
                _load_project_cached.clear()
                st.session_state["last_tree_path"]     = results["phylogeny"].get("tree_path")
                st.session_state["last_clusters"]       = results["phylogeny"].get("clusters", {})
                st.session_state["last_cluster_report"] = results["phylogeny"].get("cluster_report", {})
                st.session_state["last_matrix_path"]    = results["phylogeny"].get("matrix_path")
                st.session_state["last_upload_names"]   = results["phylogeny"].get("upload_names", [])
                st.session_state["last_cgmlst_result"]  = results["phylogeny"].get("cgmlst_result")
                if results["phylogeny"].get("error"):
                    st.warning(f"Pipeline completed with phylogeny warning: {results['phylogeny']['error']}")
                else:
                    st.success("Pipeline completed successfully.")

                if db:
                    try:
                        from backend.alerts import detect_alerts as _detect_alerts
                        _alert_data = _load_project_cached(active_project)
                        _new_alerts = _detect_alerts(_alert_data["samples"])
                        db.save_alerts(active_project, _new_alerts)
                    except Exception:
                        pass

                for _k in ("_fp_pending_job", "_fp_pending_sids", "_fp_pending_started"):
                    st.session_state.pop(_k, None)
                st.session_state["_nav_to"] = "Full Pipeline Results"
                st.rerun()

            elif _fp_job["status"] == "error":
                st.error(f"Pipeline error: {_fp_job.get('error_msg', 'Unknown error')}")
                for _k in ("_fp_pending_job", "_fp_pending_sids", "_fp_pending_started"):
                    st.session_state.pop(_k, None)

    st.divider()

    m1, m2, m3, m4 = st.columns(4)
    for col, title, tools in [
        (m1, "QC & Cleaning",    "fastp · Kraken2"),
        (m2, "QC & Assembly",  "fastp · SPAdes · Biopython"),
        (m3, "AMR Profiling",    "Minimap2 · European 2020"),
        (m4, "Phylogenetics",    "SKA2 · RapidNJ · MLST · NG-STAR"),
    ]:
        col.markdown(f"""
        <div style="
          background: #f1f5f9;
          border: 1px solid #e2e8f0;
          border-radius: 12px;
          padding: 1.25rem 1rem;
          text-align: center;
          box-shadow: 0 1px 3px rgba(0,0,0,0.06);
          min-height: 80px;
          display: flex;
          flex-direction: column;
          justify-content: center;
          gap: 0.2rem;
        ">
          <div style="font-weight:700; font-size:0.88rem; color:#64748b;">{title}</div>
          <div style="color:#64748b; font-size:0.75rem; line-height:1.5;">{tools}</div>
        </div>
        """, unsafe_allow_html=True)

    st.components.v1.html("""
<!DOCTYPE html>
<html>
<head>
<style>
  * { box-sizing:border-box; margin:0; padding:0; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif; }
  body { background:#ffffff; padding:0.6rem 0 0.4rem 0; }
  .stats-row { display:flex; gap:0; }
  .stat-item {
    flex:1; text-align:center; padding:1rem 0.5rem;
    border-right:1px solid rgba(0,0,0,0.07);
    opacity:0; transform:translateY(10px);
    transition:opacity 0.5s ease, transform 0.5s ease;
  }
  .stat-item:last-child { border-right:none; }
  .stat-item.visible { opacity:1; transform:translateY(0); }
  .stat-value {
    font-size:1.9rem; font-weight:800; color:#475569;
    line-height:1; margin-bottom:0.3rem;
    font-variant-numeric:tabular-nums; letter-spacing:-0.02em;
  }
  .stat-label {
    font-size:0.7rem; font-weight:700; color:#64748b;
    text-transform:uppercase; letter-spacing:0.09em;
  }
  .stat-sub { font-size:0.62rem; color:#64748b; margin-top:0.12rem; }
</style>
</head>
<body>
<div class="stats-row">
  <div class="stat-item" id="si0">
    <div class="stat-value" id="sv0">0</div>
    <div class="stat-label">Resistance Loci</div>
    <div class="stat-sub">20 chromosomal · 10 plasmid</div>
  </div>
  <div class="stat-item" id="si1">
    <div class="stat-value" style="font-size:1.4rem;letter-spacing:0.01em;">Genomic</div>
    <div class="stat-label">Clinical Interpretation</div>
    <div class="stat-sub">European 2020 genotype → phenotype</div>
  </div>
  <div class="stat-item" id="si2">
    <div class="stat-value">Fast</div>
    <div class="stat-label">End-to-end Runtime</div>
    <div class="stat-sub">automated · reproducible</div>
  </div>
  <div class="stat-item" id="si3">
    <div class="stat-value">CDC&nbsp;2024</div>
    <div class="stat-label">Treatment Guidelines</div>
    <div class="stat-sub">clinical interpretation</div>
  </div>
</div>
<script>
  function animCount(el, target, dur) {
    let v = 0, step = target / (dur / 16);
    function tick() {
      v = Math.min(v + step, target);
      el.textContent = Math.round(v);
      if (v < target) requestAnimationFrame(tick);
    }
    requestAnimationFrame(tick);
  }
  setTimeout(() => {
    [0,1,2,3].forEach((i, di) => {
      setTimeout(() => document.getElementById('si'+i).classList.add('visible'), di*90);
    });
    setTimeout(() => animCount(document.getElementById('sv0'), 30, 900), 120);
  }, 150);
</script>
</body>
</html>
""", height=115, scrolling=False)

    _fp_results = st.session_state.get("full_pipeline_results")
    if _fp_results:
        st.markdown("""
<div style="margin:0.5rem 0 1rem 0;padding-bottom:0.6rem;border-bottom:2px solid #e2e8f0;">
  <div style="font-size:1.05rem;font-weight:700;color:#0f172a;">Pipeline Results</div>
</div>
""", unsafe_allow_html=True)
        _all_sids = sorted(set(
            list(_fp_results.get("qc", {}).keys())
            + list(_fp_results.get("assembly", {}).keys())
            + list(_fp_results.get("amr", {}).keys())
        ))
        _fp_any_reads_qc = any(
            not _fp_results.get("qc", {}).get(_sid, {}).get("error")
            and _fp_results.get("qc", {}).get(_sid, {}).get("metrics", {}).get("total_reads")
            for _sid in _all_sids
        )

        _fp_rows, _fp_contig_paths = [], {}
        for _sid in _all_sids:
            _qr = _fp_results.get("qc",  {}).get(_sid, {})
            _ar = _fp_results.get("assembly", {}).get(_sid, {})
            _mr = _fp_results.get("amr", {}).get(_sid, {})

            _qc_reads = "—"
            _qc_gc    = "—"
            if not _qr.get("error"):
                _m = _qr.get("metrics", {})
                _tot = _m.get("total_reads")
                _gc  = _m.get("gc")
                if _tot: _qc_reads = f"{_tot:,}"
                if _gc:  _qc_gc    = f"{_gc}%"

            _asm_len = _asm_n50 = _asm_comp = "—"
            if not _ar.get("error") and _ar.get("contigs_path"):
                _st = parse_contigs_stats(_ar["contigs_path"])
                if _st is not None:
                    _asm_len  = f"{_st.total_len/1e6:.2f} Mb"
                    _asm_n50  = f"{_st.n50/1e3:.1f} kb"
                    _asm_comp = (
                        f"{_st.completeness:.1f}% ({_st.core_genes_found}/{_st.core_genes_total})"
                        if _st.core_genes_found is not None and _st.core_genes_total is not None
                        else f"{_st.completeness:.1f}%"
                    )
                    _fp_contig_paths[_sid] = Path(_ar["contigs_path"])

            _amr_score = _amr_pheno = "—"
            if isinstance(_mr, _AMRResult):
                _amr_score = f"{_mr.failure_probability*100:.1f}%"
                _amr_pheno = "; ".join(
                    c.replace("_", " ").capitalize() for c in _mr.cdc_phenotypes if c != "wildtype"
                ) or "Wildtype"

            _fp_row = {"Sample": _sid}
            if _fp_any_reads_qc:
                _fp_row["QC reads"] = _qc_reads
                _fp_row["GC%"]      = _qc_gc
            _fp_row.update({
                "Assembly length": _asm_len,
                "N50":             _asm_n50,
                "Core genes %":    _asm_comp,
                "AMR score":       _amr_score,
                "CDC phenotype":   _amr_pheno,
            })
            _fp_rows.append(_fp_row)

        _df_fp = pd.DataFrame(_fp_rows)
        st.dataframe(_df_fp, use_container_width=True, hide_index=True)

        _fpc1, _ = st.columns([2, 5])
        _fpc1.download_button(
            "⬇ Summary (CSV)",
            data=_df_fp.to_csv(index=False).encode(),
            file_name="pipeline_summary.csv",
            mime="text/csv",
            key="fp_csv_dl",
        )

        _fp_dl_cols_items = list(_fp_contig_paths.items())
        if _fp_dl_cols_items:
            st.markdown("**Download assemblies:**")
            _fpdcols = st.columns(min(len(_fp_dl_cols_items), 4))
            for _i, (_sid, _path) in enumerate(_fp_dl_cols_items):
                _fpdcols[_i % 4].download_button(
                    f"⬇ {_sid}",
                    data=_path.read_bytes(),
                    file_name=f"{_sid}_contigs.fasta",
                    mime="application/octet-stream",
                    key=f"fp_asm_dl_{_sid}",
                )

        if st.button("View full results →", type="primary", key="fp_goto_results"):
            st.session_state["_nav_to"] = "Full Pipeline Results"
            st.rerun()
        st.divider()

    st.markdown("""
<div style="margin:0.5rem 0 1rem 0;padding-bottom:0.6rem;border-bottom:2px solid #e2e8f0;">
  <div style="font-size:1.05rem;font-weight:700;color:#0f172a;">Alerts</div>
</div>
""", unsafe_allow_html=True)
    if active_project and db:
        try:
            from backend.alerts import detect_alerts as _detect_alerts_hp
            _hp_snap = _load_project_cached(active_project)
            if _hp_snap.get("samples"):
                _hp_fresh = _detect_alerts_hp(_hp_snap["samples"])
                if _hp_fresh:
                    db.save_alerts(active_project, _hp_fresh)
            _hp_all = db.load_alerts(active_project)

            _qc_errors = [a for a in _hp_all if a.get("alert_type") in {"species_contamination", "assembly_qc_fail"}]
            _mdr_xdr   = [a for a in _hp_all if a.get("alert_type") in {"mdr", "xdr"}]

            if not _hp_all:
                st.success("No alerts.")
            else:
                if _qc_errors:
                    _n_crit = sum(1 for a in _qc_errors if a.get("severity") == "critical")
                    _n_err  = sum(1 for a in _qc_errors if a.get("severity") == "error")
                    _n_warn = sum(1 for a in _qc_errors if a.get("severity") == "warning")
                    _parts  = []
                    if _n_crit: _parts.append(f"{_n_crit} critical")
                    if _n_err:  _parts.append(f"{_n_err} errors")
                    if _n_warn: _parts.append(f"{_n_warn} warnings")
                    with st.expander(f"⚠️ {' · '.join(_parts)}", expanded=False):
                        for _a in _qc_errors:
                            _sev = _a.get("severity", "")
                            _ico = "🔴" if _sev == "critical" else ("🟠" if _sev == "error" else "🟡")
                            st.markdown(f"{_ico} {_a.get('message', '')}")
                else:
                    st.success("No contamination or QC failures.")

                if _mdr_xdr:
                    _n_xdr = sum(1 for a in _mdr_xdr if a.get("alert_type") == "xdr")
                    _n_mdr = len(_mdr_xdr) - _n_xdr
                    _label = []
                    if _n_mdr:
                        _label.append(f"{_n_mdr} MDR")
                    if _n_xdr:
                        _label.append(f"{_n_xdr} XDR")
                    with st.expander(f"🔴 {' · '.join(_label)}", expanded=False):
                        for _a in sorted(_mdr_xdr, key=lambda a: a.get("alert_type") == "xdr"):
                            _icon = "☣️" if _a.get("alert_type") == "xdr" else "🔴"
                            st.markdown(f"{_icon} {_a.get('message', '')}")
                else:
                    st.success("No MDR/XDR profiles.")

        except Exception as _hp_ex:
            st.info(f"Alert system unavailable: {_hp_ex}")
    else:
        st.info("Select a project to view alerts.")

    st.divider()

    st.markdown("""
<div style="margin:0.5rem 0 0.4rem 0;padding-bottom:0.6rem;border-bottom:2px solid #e2e8f0;">
  <div style="font-size:1.05rem;font-weight:700;color:#0f172a;">Gonorrhoea Epidemiology in Europe (EU/EEA)</div>
</div>
""", unsafe_allow_html=True)

    try:
        import plotly.express as px
        import plotly.graph_objects as go
        from backend.ecdc_data import (
            EU_EEA_TOTALS   as _EU_TOTALS,
            YEARS           as _ECDC_YEARS,
            get_map_data    as _ecdc_map_data,
            SOURCE          as _ECDC_SRC,
        )

        _ISO2_ISO3 = {
            "AT": "AUT", "BE": "BEL", "BG": "BGR", "HR": "HRV", "CY": "CYP",
            "CZ": "CZE", "DK": "DNK", "EE": "EST", "FI": "FIN", "FR": "FRA",
            "DE": "DEU", "GR": "GRC", "HU": "HUN", "IS": "ISL", "IE": "IRL",
            "IT": "ITA", "LV": "LVA", "LI": "LIE", "LT": "LTU", "LU": "LUX",
            "MT": "MLT", "NL": "NLD", "NO": "NOR", "PL": "POL", "PT": "PRT",
            "RO": "ROU", "SK": "SVK", "SI": "SVN", "ES": "ESP", "SE": "SWE",
        }

        _map_col, _trend_col = st.columns([3, 1], gap="medium")

        with _map_col:
            _sel_year = st.select_slider(
                "Year",
                options=_ECDC_YEARS,
                value=2024,
                label_visibility="collapsed",
            )

            _yr_total = _EU_TOTALS[_sel_year]
            _prev     = _EU_TOTALS.get(_sel_year - 1)
            _m1, _m2, _m3 = st.columns(3)
            _m1.metric("EU/EEA Cases", f"{_yr_total['cases']:,}")
            _m2.metric("EU/EEA Rate (per 100k)", f"{_yr_total['rate']:.1f}")
            if _prev:
                _delta_r = _yr_total["rate"] - _prev["rate"]
                _m3.metric("vs Prior Year", f"{_delta_r:+.1f}")

            _raw = _ecdc_map_data(_sel_year)
            _rows = []
            for _d in _raw:
                _iso3 = _ISO2_ISO3.get(_d["iso2"])
                if not _iso3:
                    continue
                if _d["no_data"]:
                    _surv = "No data reported"
                elif _d["sentinel"]:
                    _surv = "Sentinel (cases only — no population rate)"
                else:
                    _surv = "Full national surveillance"
                _rows.append({
                    "iso_alpha":     _iso3,
                    "Country":       _d["name"],
                    "rate":          _d["rate"],
                    "Cases":         f"{_d['cases']:,}" if _d["cases"] is not None else "—",
                    "Rate per 100k": f"{_d['rate']:.1f}" if _d["rate"] is not None else "N/A",
                    "Surveillance":  _surv,
                })
            _df_map = pd.DataFrame(_rows)
            _max_r  = _df_map["rate"].max() or 120

            _fig_map = px.choropleth(
                _df_map,
                locations="iso_alpha",
                color="rate",
                hover_name="Country",
                hover_data={
                    "iso_alpha":     False,
                    "rate":          False,
                    "Cases":         True,
                    "Rate per 100k": True,
                    "Surveillance":  True,
                },
                color_continuous_scale="Reds",
                range_color=[0, max(_max_r, 120)],
                scope="europe",
                labels={"rate": "Cases per 100k"},
            )
            _fig_map.update_layout(
                paper_bgcolor="#ffffff",
                margin={"r": 0, "t": 8, "l": 0, "b": 0},
                coloraxis_colorbar=dict(
                    title="Per 100k",
                    thickness=13,
                    len=0.55,
                    tickvals=[0, 25, 50, 75, 100, 125],
                ),
                geo=dict(
                    showframe=False,
                    showcoastlines=True,
                    coastlinecolor="#d1d5db",
                    showland=True,
                    landcolor="#f3f4f6",
                    showocean=True,
                    oceancolor="#e0f2fe",
                    showlakes=False,
                    lataxis_range=[34, 72],
                    lonaxis_range=[-25, 45],
                ),
                height=440,
            )
            st.plotly_chart(_fig_map, use_container_width=True)
            st.caption(
                f"Notification rate per 100,000 population · {_ECDC_SRC['report']} · "
                f"{_ECDC_SRC['publisher']}. "
                "Grey = no rate calculated: Austria and Germany did not report; "
                "Belgium, France and Netherlands use sentinel surveillance (cases reported, "
                "no population denominator). "
                f"Data retrieved {_ECDC_SRC['retrieved']}."
            )

        with _trend_col:
            st.markdown(
                '<div style="font-size:0.7rem;font-weight:700;text-transform:uppercase;'
                'letter-spacing:0.1em;color:#64748b;margin-bottom:0.6rem;margin-top:0.4rem;">'
                'EU/EEA Trend</div>',
                unsafe_allow_html=True,
            )
            _trend_df = pd.DataFrame([
                {"Year": y, "Rate": _EU_TOTALS[y]["rate"], "Cases": _EU_TOTALS[y]["cases"]}
                for y in _ECDC_YEARS
            ])
            _fig_trend = go.Figure()
            _fig_trend.add_trace(go.Scatter(
                x=_trend_df["Year"],
                y=_trend_df["Rate"],
                mode="lines+markers",
                line=dict(color="#00c9b1", width=2),
                marker=dict(size=6, color="#00c9b1"),
                hovertemplate="%{x}: %{y:.1f} per 100k<extra></extra>",
            ))
            _fig_trend.update_layout(
                margin=dict(l=0, r=0, t=10, b=0),
                height=180,
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)",
                xaxis=dict(showgrid=False, tickfont=dict(size=10, color="#94a3b8"),
                           tickmode="array", tickvals=_ECDC_YEARS),
                yaxis=dict(showgrid=True, gridcolor="rgba(0,0,0,0.07)",
                           tickfont=dict(size=10, color="#94a3b8"), title=""),
                showlegend=False,
            )
            st.plotly_chart(_fig_trend, use_container_width=True)

            st.markdown(
                '<div style="font-size:0.7rem;font-weight:700;text-transform:uppercase;'
                'letter-spacing:0.1em;color:#64748b;margin-bottom:0.5rem;margin-top:0.6rem;">'
                f'Top 5 · {_sel_year}</div>',
                unsafe_allow_html=True,
            )
            _ranked = sorted(
                [r for r in _rows if r["rate"] is not None],
                key=lambda x: x["rate"],
                reverse=True,
            )[:5]
            for _rank, _r in enumerate(_ranked, 1):
                st.markdown(
                    f'<div style="display:flex;justify-content:space-between;'
                    f'font-size:0.75rem;margin-bottom:0.25rem;">'
                    f'<span style="color:#334155;">{_rank}. {_r["Country"]}</span>'
                    f'<span style="color:#64748b;font-weight:600;">{_r["Rate per 100k"]}</span>'
                    f'</div>',
                    unsafe_allow_html=True,
                )

    except ImportError as _imp_err:
        st.warning(f"ecdc_data module not found: {_imp_err}")
    except Exception as _map_err:
        st.info(f"Map unavailable: {_map_err}")


