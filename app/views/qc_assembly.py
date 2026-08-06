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
    _load_project_cached, _get_fastqc_html_paths,
    _inline_metadata_widget, _render_genome_overview,
    _parse_contigs_stats_backend,
)


def render() -> None:
    st.html("""
<div style="padding:0.6rem 0 1rem 0;">
  <div style="display:inline-flex;align-items:center;background:rgba(0,0,0,0.04);border:1px solid rgba(0,0,0,0.12);
  border-radius:20px;padding:0.22rem 0.85rem;font-size:0.7rem;font-weight:700;color:#64748b;
  letter-spacing:0.12em;text-transform:uppercase;margin-bottom:0.6rem;">Module 1</div>
  <div style="font-size:1.6rem;font-weight:800;color:#0f172a;">1. Quality Control &amp; Assembly</div>
</div>
""")

    _qc_proj = st.session_state.get("active_project")
    if db and _qc_proj and st.session_state.get("_qc_loaded_project") != _qc_proj:
        _qc_db = _load_project_cached(_qc_proj)
        _qc_from_db = {sid: d["qc"] for sid, d in _qc_db["samples"].items() if "qc" in d}
        st.session_state["qc_manual_results"] = _qc_from_db
        st.session_state["_qc_loaded_project"] = _qc_proj

    if _qc_proj:
        st.caption(f"Project: **{_qc_proj}** — results saved automatically")
    else:
        st.caption("No project selected — results will not be saved. Open a project on the Home page first.")

    st.markdown("""
    This module covers the full pre-assembly and assembly workflow for raw *Neisseria gonorrhoeae*
    sequencing reads — from quality control to assembled contigs ready for AMR profiling.

    **Quality Control**
    - **FastQC** — per-base quality scores, GC content, adapter content and duplication levels
    - **fastp** — adapter removal and low-quality base trimming (sliding window, minimum length)
    - **Kraken2** — species confirmation on trimmed reads; flags contaminated samples before assembly
    - **MultiQC** — aggregated HTML report summarising all QC metrics

    **De novo Assembly**
    - **SPAdes** — de Bruijn graph assembly of trimmed reads into contigs
    - **QUAST** — assembly quality metrics (N50, misassemblies) validated against *N. gonorrhoeae* thresholds
    - **chewBBACA** — core-gene completeness, % of the 1638-locus cgMLST scheme found in the assembly

    Upload raw FASTQ reads (SE or PE) to begin. After QC, select samples to proceed to assembly.
    """)

    with st.expander("Quality Control for Assembled Genomes (FASTA)"):
        st.caption(
            "Already have assembled genomes? Upload FASTA contigs here for Kraken2 species "
            "confirmation and genome statistics, no trimming or assembly needed. AMR profiling "
            "runs automatically afterwards."
        )
        _fqc_proj  = st.session_state.get("active_project")
        _fqc_files = st.file_uploader(
            "Upload assembled genomes — .fasta, .fa, .fna (one or more)",
            type=["fasta", "fa", "fna"],
            accept_multiple_files=True,
            key="fasta_qc_upload",
        )
        if _fqc_files:
            st.markdown(f"**{len(_fqc_files)} genome(s) detected**:")
            _fqc_id_map: dict = {}
            for _i, _f in enumerate(_fqc_files):
                _default_id = sanitize_name(Path(_f.name).stem)
                _c1, _c2 = st.columns([3, 2])
                _c1.markdown(f"`{_f.name}`")
                _inp = _c2.text_input(
                    "ID", value=_default_id, placeholder="e.g. NG-001",
                    key=f"fqc_sid_{_i}", label_visibility="collapsed",
                )
                _fqc_id_map[_f.name] = sanitize_name(_inp) if _inp.strip() else _default_id

            if not _fqc_proj:
                st.warning("No project selected — results will be saved to **_default**.")

            if st.button("Run QC (species confirmation + genome stats)", type="primary", key="fqc_run"):
                if db is None:
                    st.error("Database backend not available.")
                else:
                    _proj_key = _fqc_proj or "_default"
                    _fqc_job_ids: list[str] = []
                    for _f in _fqc_files:
                        _sid = _fqc_id_map[_f.name]
                        _upload_dir = PROJECTS_DIR / _proj_key / "uploads" / _sid
                        _upload_dir.mkdir(parents=True, exist_ok=True)
                        _dst = _upload_dir / _f.name
                        _dst.write_bytes(_f.getbuffer())
                        _jid = db.submit_job(
                            _proj_key, _sid, "qc_fasta",
                            {"fasta_path": str(_dst), "output_dir": str(RESULTS_BASE / "qc")},
                        )
                        _fqc_job_ids.append(_jid)
                    st.session_state["_fqc_pending_jobs"]   = _fqc_job_ids
                    st.session_state["_fqc_pending_project"] = _proj_key
                    st.session_state["fasta_qc_results"]    = {}
                    st.rerun()

        _fqc_pending = st.session_state.get("_fqc_pending_jobs", [])
        if _fqc_pending and db:
            _fqc_poll_proj = st.session_state.get("_fqc_pending_project", _fqc_proj or "")
            _all_fqc_jobs  = {j["id"]: j for j in db.get_project_jobs(_fqc_poll_proj)}
            _rel_fqc = [_all_fqc_jobs[jid] for jid in _fqc_pending if jid in _all_fqc_jobs]
            _n_fqc_tot = len(_rel_fqc)
            _n_fqc_fin = sum(1 for j in _rel_fqc if j["status"] in ("done", "error"))
            if _n_fqc_fin < _n_fqc_tot:
                st.progress(
                    _n_fqc_fin / _n_fqc_tot if _n_fqc_tot else 0,
                    text=f"{_n_fqc_fin}/{_n_fqc_tot} done",
                )
                time.sleep(2)
                st.rerun()
            else:
                _fqc_collected: dict = {}
                for _j in _rel_fqc:
                    _sid_j = _j["sample_id"]
                    if _j["status"] == "done" and _j.get("result"):
                        try:
                            _fqc_collected[_sid_j] = json.loads(_j["result"])
                        except Exception:
                            _fqc_collected[_sid_j] = {"error": "Could not parse result"}
                    else:
                        _fqc_collected[_sid_j] = {"error": _j.get("error_msg") or "Job failed"}
                st.session_state["fasta_qc_results"]  = _fqc_collected
                st.session_state["_fqc_pending_jobs"] = []
                _load_project_cached.clear()
                st.rerun()

        _fqc_results = st.session_state.get("fasta_qc_results", {})
        if _fqc_results:
            st.markdown("##### Genome QC Results")
            _fqc_rows = []
            for _sid, _res in _fqc_results.items():
                if "error" in _res:
                    _fqc_rows.append({
                        "Sample": _sid, "Species status": "ERROR", "N. gonorrhoeae": "—",
                        "Total length": "—", "Contigs": "—", "N50": "—", "AMR category": "—",
                    })
                    continue
                _k2 = (_res.get("qc") or {}).get("kraken2", {})
                _sp_status = _k2.get("species_status", "—")
                _sp_icon = {"confirmed": "✅", "likely": "⚠️", "contaminated": "❌",
                            "db_not_found": "—", "error": "⚠️"}.get(_sp_status, "—")
                _st = _res.get("assembly_stats", {})
                _amr = _res.get("amr", {})
                _fqc_rows.append({
                    "Sample":          _sid,
                    "Species status":  f"{_sp_icon} {_sp_status.replace('_', ' ').capitalize()}",
                    "N. gonorrhoeae":  f"{_k2.get('ng_pct')}%" if _k2.get("ng_pct") is not None else "—",
                    "Total length":    f"{_st['total_len'] / 1e6:.2f} Mb" if _st.get("total_len") else "—",
                    "Contigs":         _st.get("n_contigs", "—"),
                    "N50":             f"{_st['n50'] / 1e3:.1f} kb" if _st.get("n50") else "—",
                    "AMR category":    _amr.get("resistance_category", "—") if "error" not in _amr else "ERROR",
                })
            st.dataframe(pd.DataFrame(_fqc_rows), use_container_width=True, hide_index=True)
            st.download_button(
                "⬇ Summary (CSV)",
                data=pd.DataFrame(_fqc_rows).to_csv(index=False).encode(),
                file_name="fasta_qc_summary.csv",
                mime="text/csv",
                key="fqc_csv_dl",
            )

    raw_files = st.file_uploader(
        "Upload raw FASTQ files — .fastq, .fq, .fastq.gz, .fq.gz (one or more samples)",
        type=["fastq", "fq", "gz"],
        accept_multiple_files=True,
    )

    if raw_files:
        _qc_groups = group_paired_end(raw_files)
        st.markdown(f"**{len(_qc_groups)} sample(s) detected** — edit IDs if needed:")
        _qc_id_map: dict = {}
        for _i, _base in enumerate(_qc_groups):
            _is_pe = bool(_qc_groups[_base].get("R1") and _qc_groups[_base].get("R2"))
            _mode_badge = "**PE**" if _is_pe else "**SE**"
            _f_lbl = f"{_mode_badge} — " + "  +  ".join(
                f"`{_qc_groups[_base][r].name}`"
                for r in ("R1", "R2") if _qc_groups[_base].get(r)
            )
            _c1, _c2 = st.columns([3, 2])
            _c1.markdown(_f_lbl)
            _inp = _c2.text_input(
                "ID", value=_base, placeholder="e.g. NG-001",
                key=f"qc_sid_{_i}", label_visibility="collapsed",
            )
            _qc_id_map[_base] = sanitize_name(_inp) if _inp.strip() else _base

        _inline_metadata_widget(
            list(_qc_id_map.values()), _qc_proj, "qc"
        )

        _n_qc = len(_qc_groups)
        _workers_qc = min(_n_qc, max(1, os.cpu_count() or 4))
        st.caption(f"{_n_qc} sample(s) · {_workers_qc} parallel worker(s)")

        if not _qc_proj:
            st.warning("No project selected — results will be saved to **_default**. Select a project in the sidebar to save to a named project.")

        if st.button("Run QC Pipeline", type="primary"):
            if db is None:
                st.error("Database backend not available.")
            else:
                _job_ids_qc: list[str] = []
                _proj_key = _qc_proj or "_default"
                for _base, _pair in _qc_groups.items():
                    _sid = _qc_id_map[_base]
                    _upload_dir = PROJECTS_DIR / _proj_key / "uploads" / _sid
                    _upload_dir.mkdir(parents=True, exist_ok=True)
                    _file_paths: list[str] = []
                    for _rk in ("R1", "R2"):
                        if _pair.get(_rk):
                            _uf = _pair[_rk]
                            _dst = _upload_dir / _uf.name
                            _dst.write_bytes(_uf.read())
                            _file_paths.append(str(_dst))
                    _jid = db.submit_job(
                        _proj_key, _sid, "qc",
                        {"files": _file_paths, "output_dir": str(RESULTS_BASE / "qc")},
                    )
                    _job_ids_qc.append(_jid)
                st.session_state["_qc_pending_jobs"] = _job_ids_qc
                st.session_state["_qc_pending_project"] = _proj_key
                st.session_state["qc_manual_results"] = {}
                st.rerun()

    # QC job queue polling
    _qc_pending_ids = st.session_state.get("_qc_pending_jobs", [])
    if _qc_pending_ids and db:
        _qc_poll_proj = st.session_state.get("_qc_pending_project", _qc_proj or "")
        _all_qc_jobs = {j["id"]: j for j in db.get_project_jobs(_qc_poll_proj)}
        _rel_qc = [_all_qc_jobs[jid] for jid in _qc_pending_ids if jid in _all_qc_jobs]
        _n_qc_tot = len(_rel_qc)
        _n_qc_fin = sum(1 for j in _rel_qc if j["status"] in ("done", "error"))
        _n_qc_run = sum(1 for j in _rel_qc if j["status"] == "running")
        _n_qc_q   = sum(1 for j in _rel_qc if j["status"] == "queued")

        if _n_qc_fin < _n_qc_tot:
            st.markdown("#### QC running…")
            st.progress(
                _n_qc_fin / _n_qc_tot if _n_qc_tot else 0,
                text=f"{_n_qc_fin}/{_n_qc_tot} done · {_n_qc_run} running · {_n_qc_q} queued",
            )
            _running_sids = [_j["sample_id"] for _j in _rel_qc if _j["status"] == "running"]
            _queued_sids  = [_j["sample_id"] for _j in _rel_qc if _j["status"] == "queued"]
            _done_sids    = [_j["sample_id"] for _j in _rel_qc if _j["status"] in ("done", "error")]
            if _done_sids:
                st.caption(f"Done: {', '.join(_done_sids)}")
            if _running_sids:
                _dots = "." * ((int(time.time()) % 3) + 1)
                st.caption(f"Running {', '.join(_running_sids)}{_dots}")
            if _queued_sids:
                st.caption(f"Queued: {', '.join(_queued_sids)}")
            time.sleep(2)
            st.rerun()
        else:
            _qc_collected: dict = {}
            for _j in _rel_qc:
                _sid_j = _j["sample_id"]
                if _j["status"] == "done" and _j.get("result"):
                    try:
                        _qc_collected[_sid_j] = json.loads(_j["result"])
                    except Exception:
                        _qc_collected[_sid_j] = {"error": "Could not parse result"}
                else:
                    _qc_collected[_sid_j] = {"error": _j.get("error_msg") or "Job failed"}
            st.session_state["qc_manual_results"] = _qc_collected
            st.session_state["_qc_pending_jobs"] = []
            st.session_state["_qc_loaded_project"] = _qc_poll_proj
            _load_project_cached.clear()
            st.rerun()

    _qc_results = st.session_state.get("qc_manual_results", {})
    if _qc_results:
        st.markdown("---")
        st.markdown("#### QC Results")
        _qc_rows, _mqc_files = [], {}
        for _sid, _res in _qc_results.items():
            if "error" in _res:
                _qc_rows.append({
                    "Sample": _sid, "Total reads": "ERROR", "GC content": "—",
                    "N. gonorrhoeae": "—", "Species status": "—",
                })
            else:
                _met = _res.get("metrics", {})
                _k2  = _res.get("kraken2", {})
                _tot, _gc = _met.get("total_reads"), _met.get("gc")
                _ng_pct   = _k2.get("ng_pct")
                _sp_status = _k2.get("species_status", "—")
                _sp_icon = {"confirmed": "✅", "likely": "⚠️", "contaminated": "❌",
                            "db_not_found": "—", "error": "⚠️"}.get(_sp_status, "—")
                _mqp = Path(_res.get("multiqc", "") or "")
                _qc_rows.append({
                    "Sample":         _sid,
                    "Total reads":    f"{_tot:,}" if _tot else "—",
                    "GC content":     f"{_gc}%" if _gc else "—",
                    "N. gonorrhoeae": f"{_ng_pct}%" if _ng_pct is not None else "—",
                    "Species status": f"{_sp_icon} {_sp_status.replace('_', ' ').capitalize()}",
                })
                if _mqp.is_file():
                    _mqc_files[_sid] = _mqp

        _df_qc = pd.DataFrame(_qc_rows)
        st.dataframe(_df_qc, use_container_width=True, hide_index=True)

        #  Sequencing QC thresholds 
        _QC_THRESHOLDS = {
            "mean_coverage":   ("Mean coverage",    30.0,  "×",  "≥30× recommended for reliable assembly"),
            "pct_breadth_10x": ("Breadth at ≥10×",  80.0,  "%",  "≥80% of genome covered at ≥10×"),
        }
        _thresh_rows = []
        for _sid, _res in _qc_results.items():
            if "error" in _res:
                continue
            _met = _res.get("metrics", {})
            for _key, (_label, _thr, _unit, _desc) in _QC_THRESHOLDS.items():
                _val = _met.get(_key)
                if _val is None:
                    continue
                _pass = _val >= _thr
                _thresh_rows.append({
                    "Sample":    _sid,
                    "Metric":    _label,
                    "Value":     f"{_val:.1f}{_unit}",
                    "Threshold": f"≥{_thr:.0f}{_unit}",
                    "Status":    "Pass" if _pass else "Fail",
                    "Note":      _desc,
                    "_pass":     _pass,
                })

        if _thresh_rows:
            st.markdown("##### Sequencing QC — Thresholds")
            _df_thr = pd.DataFrame(_thresh_rows).drop(columns=["_pass"])

            def _thr_style(row):
                if row["Status"] == "Pass":
                    return [""] * 4 + ["background-color:#f0fdf4;color:#166534;font-weight:600"] + [""]
                return [""] * 4 + ["background-color:#fef2f2;color:#991b1b;font-weight:600"] + [""]

            st.dataframe(
                _df_thr.style.apply(_thr_style, axis=1),
                use_container_width=True,
                hide_index=True,
            )

        st.download_button(
            "⬇ Download summary (CSV)",
            data=_df_qc.to_csv(index=False).encode(),
            file_name="qc_summary.csv",
            mime="text/csv",
            key="qc_csv_dl",
        )

        for _sid, _res in _qc_results.items():
            if "error" in _res:
                continue
            _fqc_paths = _get_fastqc_html_paths(_res)
            for _fqc_html in _fqc_paths:
                _fqc_label = _fqc_html.stem.replace("_fastqc", "")
                with st.expander(f"FastQC Report — {_fqc_label}", expanded=False):
                    _fqc_col1, _fqc_col2 = st.columns([5, 1])
                    with _fqc_col1:
                        st.components.v1.html(
                            _fqc_html.read_text(encoding="utf-8"),
                            height=800,
                            scrolling=True,
                        )
                    with _fqc_col2:
                        st.download_button(
                            "⬇ Download HTML",
                            data=_fqc_html.read_bytes(),
                            file_name=_fqc_html.name,
                            mime="text/html",
                            key=f"dl_fastqc_{_fqc_html.stem}",
                        )
            _mqc_path = _res.get("multiqc")
            if _mqc_path and Path(_mqc_path).exists():
                with st.expander(f"MultiQC Summary — {_sid}", expanded=False):
                    _mqc_col1, _mqc_col2 = st.columns([5, 1])
                    with _mqc_col1:
                        st.components.v1.html(
                            Path(_mqc_path).read_text(encoding="utf-8"),
                            height=900,
                            scrolling=True,
                        )
                    with _mqc_col2:
                        st.download_button(
                            "⬇ Download MultiQC",
                            data=Path(_mqc_path).read_bytes(),
                            file_name="multiqc_report.html",
                            mime="text/html",
                            key=f"dl_multiqc_{_sid}",
                        )

        # Checkpoint: proceed to assembly
        _samples_with_trimmed = {
            sid: res for sid, res in _qc_results.items()
            if "error" not in res and res.get("trimmed_r1")
        }
        if _samples_with_trimmed:
            st.markdown("---")
            st.markdown("#### Proceed to Assembly")
            st.caption("Samples with QC issues are pre-deselected. Untick any sample you want to exclude before running assembly.")

            def _qc_verdict(res):
                sp = res.get("kraken2", {}).get("species_status", "")
                reads = res.get("metrics", {}).get("total_reads") or 0
                if sp == "contaminated":
                    return "fail", "Contaminated"
                if reads < 10_000:
                    return "fail", f"Too few reads ({reads:,})"
                if sp == "likely" or reads < 50_000:
                    return "warn", "Low quality — review"
                return "pass", "Pass"

            _asm_selected: dict[str, bool] = {}
            for _sid, _res in _samples_with_trimmed.items():
                _verdict, _reason = _qc_verdict(_res)
                _default_sel = _verdict != "fail"
                _asm_selected[_sid] = st.checkbox(
                    _sid, value=_default_sel, key=f"qc2asm_{_sid}"
                )

            _n_sel = sum(_asm_selected.values())
            if st.button(
                f"Run Assembly — {_n_sel} sample(s) selected",
                type="primary",
                disabled=_n_sel == 0,
                key="qc2asm_run",
            ):
                if db is None:
                    st.error("Database backend not available.")
                else:
                    _to_asm = {
                        sid: res for sid, res in _samples_with_trimmed.items()
                        if _asm_selected.get(sid)
                    }
                    _asm_job_ids: list[str] = []
                    _proj_key_asm = _qc_proj or "_default"
                    for _sid, _res in _to_asm.items():
                        _jid = db.submit_job(
                            _proj_key_asm, _sid, "assembly",
                            {
                                "r1": _res["trimmed_r1"],
                                "r2": _res.get("trimmed_r2"),
                                "pre_trimmed": True,
                            },
                        )
                        _asm_job_ids.append(_jid)
                    st.session_state["_asm_pending_jobs"] = _asm_job_ids
                    st.session_state["_asm_pending_project"] = _proj_key_asm
                    st.session_state["assembly_results"] = {}
                    st.rerun()

        # Restore saved assembly results from DB
        _asm_proj = st.session_state.get("active_project")
        if db and _asm_proj and st.session_state.get("_asm_loaded_project") != _asm_proj:
            _asm_db = _load_project_cached(_asm_proj)
            _asm_from_db = {
                sid: {"contigs_path": d["assembly"]["contigs_path"]}
                for sid, d in _asm_db["samples"].items()
                if "assembly" in d and d["assembly"].get("contigs_path")
            }
            st.session_state["assembly_results"] = _asm_from_db
            st.session_state["_asm_loaded_project"] = _asm_proj

    # Assembly job queue polling
    _asm_pending_ids = st.session_state.get("_asm_pending_jobs", [])
    if _asm_pending_ids and db:
        _asm_poll_proj = st.session_state.get("_asm_pending_project", _qc_proj or "")
        _all_asm_jobs = {j["id"]: j for j in db.get_project_jobs(_asm_poll_proj)}
        _rel_asm = [_all_asm_jobs[jid] for jid in _asm_pending_ids if jid in _all_asm_jobs]
        _n_atot = len(_rel_asm)
        _n_afin = sum(1 for j in _rel_asm if j["status"] in ("done", "error"))
        _n_arun = sum(1 for j in _rel_asm if j["status"] == "running")
        _n_aq   = sum(1 for j in _rel_asm if j["status"] == "queued")

        if _n_afin < _n_atot:
            st.markdown("#### Assembly running…")
            st.progress(
                _n_afin / _n_atot if _n_atot else 0,
                text=f"{_n_afin}/{_n_atot} done · {_n_arun} running",
            )
            _asm_done_sids    = [_j["sample_id"] for _j in _rel_asm if _j["status"] in ("done", "error")]
            _asm_running_sids = [_j["sample_id"] for _j in _rel_asm if _j["status"] == "running"]
            if _asm_done_sids:
                st.caption(f"Done: {', '.join(_asm_done_sids)}")
            if _asm_running_sids:
                _dots = "." * ((int(time.time()) % 3) + 1)
                st.caption(f"Running {', '.join(_asm_running_sids)}{_dots}")
            time.sleep(2)
            st.rerun()
        else:
            _asm_collected: dict = {}
            _amr_job_ids: list[str] = []
            for _j in _rel_asm:
                _sid_j = _j["sample_id"]
                if _j["status"] == "done" and _j.get("result"):
                    try:
                        _asm_collected[_sid_j] = json.loads(_j["result"])
                    except Exception:
                        _asm_collected[_sid_j] = {"error": "Could not parse result"}
                else:
                    _asm_collected[_sid_j] = {"error": _j.get("error_msg") or "Job failed"}

            # Auto-submit AMR for every successfully assembled sample
            for _sid_j, _ares in _asm_collected.items():
                _cpath = _ares.get("contigs_path")
                if _cpath and "error" not in _ares and db:
                    _amr_payload: dict = {"contigs": _cpath}
                    _amr_job_ids.append(
                        db.submit_job(_asm_poll_proj, _sid_j, "amr", _amr_payload)
                    )

            st.session_state["assembly_results"] = _asm_collected
            st.session_state["_asm_pending_jobs"] = []
            if _amr_job_ids:
                st.session_state["_amr_worker_pending_jobs"] = _amr_job_ids
                st.session_state["_amr_worker_pending_project"] = _asm_poll_proj
            _load_project_cached.clear()
            st.rerun()

    # AMR worker job polling (auto-triggered after assembly)
    _amr_worker_ids = st.session_state.get("_amr_worker_pending_jobs", [])
    if _amr_worker_ids and db:
        _amr_wpoll_proj = st.session_state.get("_amr_worker_pending_project", "")
        _all_amrw_jobs = {j["id"]: j for j in db.get_project_jobs(_amr_wpoll_proj)}
        _rel_amrw = [_all_amrw_jobs[jid] for jid in _amr_worker_ids if jid in _all_amrw_jobs]
        _n_amrw_tot = len(_rel_amrw)
        _n_amrw_fin = sum(1 for j in _rel_amrw if j["status"] in ("done", "error"))
        _n_amrw_run = sum(1 for j in _rel_amrw if j["status"] == "running")
        _n_amrw_q   = sum(1 for j in _rel_amrw if j["status"] == "queued")
        if _n_amrw_fin < _n_amrw_tot:
            st.markdown("#### AMR analysis running…")
            st.progress(
                _n_amrw_fin / _n_amrw_tot if _n_amrw_tot else 0,
                text=f"{_n_amrw_fin}/{_n_amrw_tot} done · {_n_amrw_run} running · {_n_amrw_q} queued",
            )
            for _j in _rel_amrw:
                _ico = {"queued": "⏳", "running": "⏩", "done": "✓", "error": "✗"}.get(_j["status"], "·")
                st.caption(f"{_ico} {_j['sample_id']} — {_j['status']}")
            time.sleep(2)
            st.rerun()
        else:
            st.session_state["_amr_worker_pending_jobs"] = []
            _load_project_cached.clear()
            st.rerun()

    # Full assembly results (from QC checkpoint or direct assembly below)
    _asm_results = st.session_state.get("assembly_results", {})

    # Direct Assembly expander (for pre-trimmed reads or skipping QC)
    with st.expander("Direct Assembly — assemble reads without running QC first"):
        st.caption("Use this when reads are already trimmed or you want to skip the QC step.")
        _da_proj = st.session_state.get("active_project")
        _da_uploaded = st.file_uploader(
            "Upload FASTQ files — .fastq, .fq, .fastq.gz (one or more samples)",
            type=["fastq", "fq", "gz"],
            accept_multiple_files=True,
            key="da_upload",
        )
        if _da_uploaded:
            _da_groups = group_paired_end(_da_uploaded)
            st.markdown(f"**{len(_da_groups)} sample(s) detected** — edit IDs if needed:")
            _da_id_map: dict = {}
            for _i, _base in enumerate(_da_groups):
                _is_pe = bool(_da_groups[_base].get("R1") and _da_groups[_base].get("R2"))
                _mode_badge = "**PE**" if _is_pe else "**SE**"
                _f_lbl = f"{_mode_badge} — " + "  +  ".join(
                    f"`{_da_groups[_base][r].name}`"
                    for r in ("R1", "R2") if _da_groups[_base].get(r)
                )
                _c1, _c2 = st.columns([3, 2])
                _c1.markdown(_f_lbl)
                _inp = _c2.text_input(
                    "ID", value=_base, placeholder="e.g. NG-001",
                    key=f"da_sid_{_i}", label_visibility="collapsed",
                )
                _da_id_map[_base] = sanitize_name(_inp) if _inp.strip() else _base

            _inline_metadata_widget(list(_da_id_map.values()), _da_proj, "da")

            if st.button("Run Assembly", key="da_run_btn"):
                if run_assembly_pipeline is None:
                    st.error("Assembly backend not available.")
                else:
                    _da_all: dict = {}
                    _n_da = len(_da_groups)
                    _workers_da = max(1, min(_n_da, (os.cpu_count() or 2) // 2))
                    _prog_da = st.progress(0, text="Saving uploads…")
                    _da_temp: dict[str, tuple] = {}
                    _da_saved: list[Path] = []
                    for _base, _pair in _da_groups.items():
                        _sid = _da_id_map[_base]
                        _tr1 = save_temp_file(_pair["R1"]) if _pair.get("R1") else None
                        _tr2 = save_temp_file(_pair["R2"]) if _pair.get("R2") else None
                        _da_temp[_sid] = (_tr1, _tr2)
                        if _tr1: _da_saved.append(_tr1)
                        if _tr2: _da_saved.append(_tr2)
                    _prog_da.progress(0, text=f"Assembling {_n_da} sample(s)…")
                    _done_da = 0
                    try:
                        with ThreadPoolExecutor(max_workers=_workers_da) as _exec_da:
                            _da_futures = {
                                _exec_da.submit(
                                    run_assembly_pipeline, _sid, _r1, _r2,
                                ): _sid
                                for _sid, (_r1, _r2) in _da_temp.items()
                            }
                            for _fut in as_completed(_da_futures):
                                _sid = _da_futures[_fut]
                                _done_da += 1
                                _prog_da.progress(_done_da / _n_da, text=f"Done {_done_da}/{_n_da} — {_sid}")
                                try:
                                    _asm = _fut.result()
                                    if _asm.contigs:
                                        _da_all[_sid] = {
                                            "contigs_path": str(_asm.contigs),
                                            "logs": _asm.logs,
                                        }
                                    else:
                                        _da_all[_sid] = {"error": "Assembly failed", "logs": _asm.logs}
                                except Exception as _e:
                                    _da_all[_sid] = {"error": str(_e)}
                    finally:
                        cleanup_temp_files(_da_saved)
                    _prog_da.progress(1.0, text="Done")
                    st.session_state["assembly_results"] = _da_all
                    if db and _da_proj:
                        db.init_project(_da_proj)
                        _da_amr_jobs: list[str] = []
                        for _sid, _res in _da_all.items():
                            if "error" not in _res and _res.get("contigs_path"):
                                from backend.phylogeny.cgmlst import core_genome_completeness
                                try:
                                    core_genome_completeness(
                                        _res["contigs_path"], str(Path(_res["contigs_path"]).parent)
                                    )
                                except Exception:
                                    pass
                                _stats = _parse_contigs_stats_backend(Path(_res["contigs_path"]))
                                db.save_assembly(_da_proj, _sid, _res["contigs_path"], _stats)
                                _da_amr_payload: dict = {"contigs": _res["contigs_path"]}
                                _da_amr_jobs.append(
                                    db.submit_job(_da_proj, _sid, "amr", _da_amr_payload)
                                )
                        if _da_amr_jobs:
                            st.session_state["_amr_worker_pending_jobs"] = _da_amr_jobs
                            st.session_state["_amr_worker_pending_project"] = _da_proj
                        _load_project_cached.clear()
                    st.rerun()

    _asm_results = st.session_state.get("assembly_results", {})
    if _asm_results:
        st.markdown("---")
        st.markdown("#### Assembly Results")
        _asm_rows, _contig_paths = [], {}
        for _sid, _res in _asm_results.items():
            if "error" in _res:
                _asm_rows.append({
                    "Sample": _sid,
                    "Status": "❌",
                    "Total length": "—",
                    "Contigs": "—",
                    "N50": "—",
                    "GC%": "—",
                    "Core genes %": "—",
                    "Misassemblies": "—",
                    "Duplication ratio": "—",
                    "Mismatches/100kbp": "—",
                    "Indels/100kbp": "—",
                    "QC flags": _res.get("error", "Error"),
                })
                if _res.get("logs"):
                    with st.expander(f"{_sid} — error logs"):
                        st.text(_res["logs"][-2000:])
            else:
                _cp = _res.get("contigs_path", "")
                _st = parse_contigs_stats(_cp) if _cp else None
                if _st is not None:
                    _status = {"pass": "✅", "warn": "⚠️", "fail": "❌"}.get(_st.qc_status, "—")
                    _flags  = "; ".join(_st.qc_flags) or "—"
                    if _st.core_genes_found is not None and _st.core_genes_total is not None:
                        _core_str = f"{_st.completeness:.1f}% ({_st.core_genes_found}/{_st.core_genes_total})"
                    else:
                        _core_str = f"{_st.completeness:.1f}%"
                    _asm_rows.append({
                        "Sample":          _sid,
                        "Status":          _status,
                        "Total length":    f"{_st.total_len / 1e6:.2f} Mb",
                        "Contigs":         _st.n_contigs,
                        "N50":             f"{_st.n50 / 1e3:.1f} kb",
                        "GC%":             f"{_st.gc_pct:.1f}%",
                        "Core genes %":    _core_str,
                        "Misassemblies":   str(_st.misassemblies) if _st.misassemblies is not None else "—",
                        "Duplication ratio": f"{_st.duplication_ratio:.2f}" if _st.duplication_ratio is not None else "—",
                        "Mismatches/100kbp": f"{_st.mismatches_per_100kbp:.2f}" if _st.mismatches_per_100kbp is not None else "—",
                        "Indels/100kbp":     f"{_st.indels_per_100kbp:.2f}" if _st.indels_per_100kbp is not None else "—",
                        "QC flags":        _flags,
                    })
                    _contig_paths[_sid] = Path(_cp)

        _df_asm = pd.DataFrame(_asm_rows)
        if not _df_asm.empty:
            def _asm_row_style(row):
                if row.get("Status") == "❌":
                    return ["background-color:#fef2f2;color:#dc2626"] * len(row)
                if row.get("Status") == "⚠️":
                    return ["background-color:#fffbeb;color:#92400e"] * len(row)
                return [""] * len(row)
            st.dataframe(
                _df_asm.style.apply(_asm_row_style, axis=1),
                use_container_width=True,
                hide_index=True,
            )
        elif not any("error" in _r for _r in _asm_results.values()):
            st.info("Assembly statistics unavailable — contig files not found.")

        st.download_button(
            "⬇ Summary (CSV)",
            data=_df_asm.to_csv(index=False).encode(),
            file_name="assembly_summary.csv",
            mime="text/csv",
            key="asm_csv_dl",
        )
        if _contig_paths:
            st.markdown("**Download assembled genomes:**")
            _dcols = st.columns(min(len(_contig_paths), 4))
            for _i, (_sid, _path) in enumerate(_contig_paths.items()):
                _dcols[_i % 4].download_button(
                    f"⬇ {_sid}",
                    data=_path.read_bytes(),
                    file_name=f"{_sid}_contigs.fasta",
                    mime="application/octet-stream",
                    key=f"asm_dl_{_sid}",
                )
            for _sid, _path in _contig_paths.items():
                st.markdown(f"**{_sid}**")
                _render_genome_overview(str(_path), {}, _sid)

