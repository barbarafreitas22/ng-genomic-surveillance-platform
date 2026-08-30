from __future__ import annotations
import json
import os
import time
from pathlib import Path

import pandas as pd
import streamlit as st

from views.shared import *
from views.shared import (
    _load_project_cached,
    _inline_metadata_widget,
)


def _qc_verdict(res: dict) -> tuple[str, str]:
    sp = res.get("kraken2", {}).get("species_status", "")
    reads = res.get("metrics", {}).get("total_reads") or 0
    if sp == "contaminated":
        return "fail", "Contaminated"
    if reads < 10_000:
        return "fail", f"Too few reads ({reads:,})"
    if sp == "likely" or reads < 50_000:
        return "warn", "Low quality — review"
    return "pass", "Pass"


def _job_status_counts(jobs: list[dict]) -> dict:
    return {
        "total":   len(jobs),
        "done":    sum(1 for j in jobs if j["status"] in ("done", "error")),
        "running": sum(1 for j in jobs if j["status"] == "running"),
        "queued":  sum(1 for j in jobs if j["status"] == "queued"),
    }


def _collect_job_results(jobs: list[dict]) -> dict:
    collected: dict = {}
    for j in jobs:
        sid = j["sample_id"]
        if j["status"] == "done" and j.get("result"):
            try:
                collected[sid] = json.loads(j["result"])
            except Exception:
                collected[sid] = {"error": "Could not parse result"}
        else:
            collected[sid] = {"error": j.get("error_msg") or "Job failed"}
    return collected


def _render_job_wait(counts: dict, heading: str | None, progress_text: str, captions: list[str]) -> None:
    if heading:
        st.markdown(heading)
    st.progress(counts["done"] / counts["total"] if counts["total"] else 0, text=progress_text)
    for c in captions:
        st.caption(c)
    time.sleep(2)
    st.rerun()


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
        st.caption(f"Project: **{_qc_proj}**")
        if st.button("Delete all results", key="del_qc_asm_btn"):
            from backend.db import delete_qc as _delete_qc, delete_assembly as _delete_assembly
            _delete_qc(_qc_proj)
            _delete_assembly(_qc_proj)
            for _k in ("qc_manual_results", "fasta_qc_results", "assembly_results",
                       "_qc_loaded_project", "_asm_loaded_project"):
                st.session_state.pop(_k, None)
            _load_project_cached.clear()
            st.rerun()
    else:
        st.caption("No project selected — results will not be saved. Open a project on the Home page first.")

    st.markdown("""
    This module covers the full pre-assembly and assembly workflow for *Neisseria gonorrhoeae*
    sequencing data. Two input types are accepted, and the pipeline runs different steps
    depending on which one is used.

    **FASTQ**: fastp and Kraken2 on the reads, then SPAdes assembly, contiguity
    metrics (Biopython) and core genome completeness on the resulting contigs.

    **FASTA**: Kraken2, contiguity metrics (Biopython) and core genome completeness run
    directly on the uploaded contigs.

    | Stage | Tool | What it does |
    |---|---|---|
    | Quality Control | **fastp** | Adapter removal and low-quality base trimming (sliding window, minimum length); pre-trimming read count and GC% also come from fastp's own report |
    | Quality Control | **Kraken2** | Species confirmation on a locally-constructed database; flags contaminated samples before assembly |
    | De novo Assembly | **SPAdes** | De Bruijn graph assembly of trimmed reads into contigs |
    | De Novo Assembly Evaluation | **Biopython** | Contiguity metrics (N50, N90, L50, L90, auN, GC%, contig count) computed directly from the assembly FASTA and validated against *N. gonorrhoeae* thresholds |
    | Core Genome Completeness | **BLASTN** | Against a local set of 1,713 *N. gonorrhoeae* core genes. Reports the % of genes found present and intact, and is checked independently of the contiguity metrics above |
    """)

    _fqc_proj = st.session_state.get("active_project")
    _FASTA_EXTS = (".fasta", ".fa", ".fna")

    _all_files = st.file_uploader(
        "Upload raw FASTQ reads or assembled genomes (FASTA). The type of each file is "
        "detected automatically from its extension.",
        type=["fastq", "fq", "gz", "fasta", "fa", "fna"],
        accept_multiple_files=True,
        key="qc_upload_unified",
    )
    _fqc_files   = [f for f in (_all_files or []) if f.name.lower().endswith(_FASTA_EXTS)]
    _fastq_files = [f for f in (_all_files or []) if f not in _fqc_files]

    if _fqc_files or _fastq_files:
        _qc_groups = group_paired_end(_fastq_files) if _fastq_files else {}
        _n_detected = len(_qc_groups) + len(_fqc_files)
        st.markdown(f"**{_n_detected} {plural(_n_detected, 'sample')} detected** — edit IDs if needed:")

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

        _fqc_id_map: dict = {}
        for _i, _f in enumerate(_fqc_files):
            _default_id = sanitize_name(Path(_f.name).stem)
            _c1, _c2 = st.columns([3, 2])
            _c1.markdown(f"**FASTA** — `{_f.name}`")
            _inp = _c2.text_input(
                "ID", value=_default_id, placeholder="e.g. NG-001",
                key=f"fqc_sid_{_i}", label_visibility="collapsed",
            )
            _fqc_id_map[_f.name] = sanitize_name(_inp) if _inp.strip() else _default_id

        _inline_metadata_widget(
            list(_qc_id_map.values()) + list(_fqc_id_map.values()), _qc_proj, "qc"
        )

        if _qc_groups:
            _n_qc = len(_qc_groups)
            _workers_qc = min(_n_qc, max(1, os.cpu_count() or 4))
            st.caption(f"{_n_qc} FASTQ {plural(_n_qc, 'sample')} · {_workers_qc} parallel {plural(_workers_qc, 'worker')}")

        if not _qc_proj:
            st.warning("No project selected. Select a project in the sidebar to save to a named project.")

        if st.button("Run QC Pipeline", type="primary"):
            if db is None:
                st.error("Database backend not available.")
            else:
                _proj_key = _qc_proj or "_default"

                _job_ids_qc: list[str] = []
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
                st.session_state["_fqc_pending_jobs"]    = _fqc_job_ids
                st.session_state["_fqc_pending_project"] = _proj_key
                st.session_state["fasta_qc_results"]     = {}

                st.rerun()

    _fqc_pending = st.session_state.get("_fqc_pending_jobs", [])
    if _fqc_pending and db:
        _fqc_poll_proj = st.session_state.get("_fqc_pending_project", _fqc_proj or "")
        _all_fqc_jobs  = {j["id"]: j for j in db.get_project_jobs(_fqc_poll_proj)}
        _rel_fqc = [_all_fqc_jobs[jid] for jid in _fqc_pending if jid in _all_fqc_jobs]
        _c_fqc = _job_status_counts(_rel_fqc)
        if _c_fqc["done"] < _c_fqc["total"]:
            _render_job_wait(_c_fqc, None, f"{_c_fqc['done']}/{_c_fqc['total']} done", [])
        else:
            _fqc_collected = _collect_job_results(_rel_fqc)
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
                    "Total length": "—", "Contigs": "—", "N50": "—",
                    "Avg contig length": "—", "N's/100kbp": "—", "Core genes %": "—",
                })
                continue
            _k2 = (_res.get("qc") or {}).get("kraken2", {})
            _sp_status = _k2.get("species_status", "—")
            _sp_icon = {"confirmed": "✅", "likely": "⚠️", "contaminated": "❌",
                        "db_not_found": "—", "error": "⚠️"}.get(_sp_status, "—")
            _st = _res.get("assembly_stats", {})
            _fqc_rows.append({
                "Sample":          _sid,
                "Species status":  f"{_sp_icon} {_sp_status.replace('_', ' ').capitalize()}",
                "N. gonorrhoeae":  f"{_k2.get('ng_pct')}%" if _k2.get("ng_pct") is not None else "—",
                "Total length":    f"{_st['total_len'] / 1e6:.2f} Mb" if _st.get("total_len") else "—",
                "Contigs":         _st.get("n_contigs", "—"),
                "N50":             f"{_st['n50'] / 1e3:.1f} kb" if _st.get("n50") else "—",
                "Avg contig length": f"{_st['avg_contig_len'] / 1e3:.1f} kb" if _st.get("avg_contig_len") else "—",
                "N's/100kbp":      f"{_st['n_per_100kbp']:.1f}" if _st.get("n_per_100kbp") is not None else "—",
                "Core genes %":    f"{_st['completeness']:.1f}%" if _st.get("completeness") else "—",
            })
        st.dataframe(pd.DataFrame(_fqc_rows), use_container_width=True, hide_index=True)
        st.download_button(
            "⬇ Summary (CSV)",
            data=pd.DataFrame(_fqc_rows).to_csv(index=False).encode(),
            file_name="fasta_qc_summary.csv",
            mime="text/csv",
            key="fqc_csv_dl",
        )

    # QC job queue polling
    _qc_pending_ids = st.session_state.get("_qc_pending_jobs", [])
    if _qc_pending_ids and db:
        _qc_poll_proj = st.session_state.get("_qc_pending_project", _qc_proj or "")
        _all_qc_jobs = {j["id"]: j for j in db.get_project_jobs(_qc_poll_proj)}
        _rel_qc = [_all_qc_jobs[jid] for jid in _qc_pending_ids if jid in _all_qc_jobs]
        _c_qc = _job_status_counts(_rel_qc)

        if _c_qc["done"] < _c_qc["total"]:
            _running_sids = [_j["sample_id"] for _j in _rel_qc if _j["status"] == "running"]
            _queued_sids  = [_j["sample_id"] for _j in _rel_qc if _j["status"] == "queued"]
            _done_sids    = [_j["sample_id"] for _j in _rel_qc if _j["status"] in ("done", "error")]
            _captions = []
            if _done_sids:
                _captions.append(f"Done: {', '.join(_done_sids)}")
            if _running_sids:
                _dots = "." * ((int(time.time()) % 3) + 1)
                _captions.append(f"Running {', '.join(_running_sids)}{_dots}")
            if _queued_sids:
                _captions.append(f"Queued: {', '.join(_queued_sids)}")
            _render_job_wait(
                _c_qc, "#### QC running…",
                f"{_c_qc['done']}/{_c_qc['total']} done · {_c_qc['running']} running · {_c_qc['queued']} queued",
                _captions,
            )
        else:
            _qc_collected = _collect_job_results(_rel_qc)
            st.session_state["qc_manual_results"] = _qc_collected
            st.session_state["_qc_pending_jobs"] = []
            st.session_state["_qc_loaded_project"] = _qc_poll_proj

            if db is not None:
                _auto_asm_ids: list[str] = []
                for _sid_j, _res_j in _qc_collected.items():
                    if "error" in _res_j or not _res_j.get("trimmed_r1"):
                        continue
                    if _qc_verdict(_res_j)[0] == "fail":
                        continue
                    _auto_asm_ids.append(db.submit_job(
                        _qc_poll_proj, _sid_j, "assembly",
                        {
                            "r1": _res_j["trimmed_r1"],
                            "r2": _res_j.get("trimmed_r2"),
                            "pre_trimmed": True,
                        },
                    ))
                if _auto_asm_ids:
                    st.session_state["_asm_pending_jobs"] = _auto_asm_ids
                    st.session_state["_asm_pending_project"] = _qc_poll_proj
                    st.session_state["assembly_results"] = {}

            _load_project_cached.clear()
            st.rerun()

    _qc_results = st.session_state.get("qc_manual_results", {})
    if _qc_results:
        st.markdown("---")
        st.markdown("#### QC Results")
        _qc_rows = []
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
                _qc_rows.append({
                    "Sample":         _sid,
                    "Total reads":    f"{_tot:,}" if _tot else "—",
                    "GC content":     f"{_gc}%" if _gc else "—",
                    "N. gonorrhoeae": f"{_ng_pct}%" if _ng_pct is not None else "—",
                    "Species status": f"{_sp_icon} {_sp_status.replace('_', ' ').capitalize()}",
                })

        _df_qc = pd.DataFrame(_qc_rows)
        st.dataframe(_df_qc, use_container_width=True, hide_index=True)

        #  Sequencing QC thresholds 
        _QC_THRESHOLDS = {
            "mean_coverage":   ("Mean coverage",    40.0,  "×",  "≥40× per CDC AR Lab Network EQA (Reimche et al. 2026)"),
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


        _already_asm = set(st.session_state.get("assembly_results", {}).keys())
        _asm_pending_now = bool(st.session_state.get("_asm_pending_jobs"))
        _samples_with_trimmed = {
            sid: res for sid, res in _qc_results.items()
            if "error" not in res and res.get("trimmed_r1") and sid not in _already_asm
        }
        if _samples_with_trimmed and not _asm_pending_now:
            st.markdown("---")
            st.markdown("#### Proceed to Assembly")
            st.caption("Samples with QC issues are pre-deselected. Untick any sample you want to exclude before running assembly.")

            _keep_spades = st.checkbox(
                "Keep SPAdes intermediate files", value=False,
                help="Preserves the assembly graph, per-k-mer directories, and spades.log "
                     "(~500 MB–2 GB per sample) instead of deleting them after the contigs FASTA is copied out.",
                key="qc2asm_keep_spades",
            )
            if _keep_spades:
                st.warning(
                    "SPAdes intermediate files take ~500 MB–2 GB per sample "
                    "(10–20 GB for 10 samples). Make sure there's enough disk space."
                )

            _asm_selected: dict[str, bool] = {}
            for _sid, _res in _samples_with_trimmed.items():
                _verdict, _reason = _qc_verdict(_res)
                _default_sel = _verdict != "fail"
                _asm_selected[_sid] = st.checkbox(
                    _sid, value=_default_sel, key=f"qc2asm_{_sid}"
                )

            _n_sel = sum(_asm_selected.values())
            if st.button(
                f"Run Assembly: {_n_sel} {plural(_n_sel, 'sample')} selected",
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
                                "keep_spades_output": _keep_spades,
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
                sid: {
                    "contigs_path": d["assembly"]["contigs_path"],
                    "is_reads_derived": d.get("qc", {}).get("metrics", {}).get("total_reads") is not None,
                }
                for sid, d in _asm_db["samples"].items()
                if "assembly" in d and d["assembly"].get("contigs_path")
            }
            st.session_state["assembly_results"] = _asm_from_db
            st.session_state["_asm_loaded_project"] = _asm_proj

    
    _asm_pending_ids = st.session_state.get("_asm_pending_jobs", [])
    if _asm_pending_ids and db:
        _asm_poll_proj = st.session_state.get("_asm_pending_project", _qc_proj or "")
        _all_asm_jobs = {j["id"]: j for j in db.get_project_jobs(_asm_poll_proj)}
        _rel_asm = [_all_asm_jobs[jid] for jid in _asm_pending_ids if jid in _all_asm_jobs]
        _c_asm = _job_status_counts(_rel_asm)

        if _c_asm["done"] < _c_asm["total"]:
            _asm_done_sids    = [_j["sample_id"] for _j in _rel_asm if _j["status"] in ("done", "error")]
            _asm_running_sids = [_j["sample_id"] for _j in _rel_asm if _j["status"] == "running"]
            _captions = []
            if _asm_done_sids:
                _captions.append(f"Done: {', '.join(_asm_done_sids)}")
            if _asm_running_sids:
                _dots = "." * ((int(time.time()) % 3) + 1)
                _captions.append(f"Running {', '.join(_asm_running_sids)}{_dots}")
            _render_job_wait(
                _c_asm, "#### Assembly running…",
                f"{_c_asm['done']}/{_c_asm['total']} done · {_c_asm['running']} running",
                _captions,
            )
        else:
            _asm_collected = _collect_job_results(_rel_asm)
            _amr_job_ids: list[str] = []
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
        _c_amrw = _job_status_counts(_rel_amrw)
        if _c_amrw["done"] < _c_amrw["total"]:
            _captions = []
            for _j in _rel_amrw:
                _ico = {"queued": "⏳", "running": "⏩", "done": "✓", "error": "✗"}.get(_j["status"], "·")
                _captions.append(f"{_ico} {_j['sample_id']}: {_j['status']}")
            _render_job_wait(
                _c_amrw, "#### AMR analysis running…",
                f"{_c_amrw['done']}/{_c_amrw['total']} done · {_c_amrw['running']} running · {_c_amrw['queued']} queued",
                _captions,
            )
        else:
            st.session_state["_amr_worker_pending_jobs"] = []
            _load_project_cached.clear()
            st.rerun()

    _asm_results = st.session_state.get("assembly_results", {})
    if _asm_results:
        st.markdown("---")
        _asm_rows, _contig_paths, _reads_derived_present = [], {}, False
        for _sid, _res in _asm_results.items():
            if not _res.get("is_reads_derived", True):
                _cp = _res.get("contigs_path", "")
                if _cp:
                    _contig_paths[_sid] = Path(_cp)
                continue
            _reads_derived_present = True
            if "error" in _res:
                _asm_rows.append({
                    "Sample": _sid,
                    "Status": "❌",
                    "Total length": "—",
                    "Contigs": "—",
                    "N50": "—",
                    "GC%": "—",
                    "N's/100kbp": "—",
                    "Core genes %": "—",
                    "QC flags": _res.get("error", "Error"),
                })
                if _res.get("logs"):
                    with st.expander(f"{_sid} — error logs"):
                        st.text(_res["logs"][-2000:])
            else:
                _cp = _res.get("contigs_path", "")
                _st = parse_contigs_stats(_cp) if _cp else None
                if _st is not None:
                    if _st.qc_status == "fail" and not any(
                        f.startswith("Core genes") or f.startswith("GC%") for f in _st.qc_flags
                    ):
                        _status = "❌ Fragmentation"
                    else:
                        _status = {"pass": "✅", "caution": "⚠️", "fail": "❌"}.get(_st.qc_status, "—")
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
                        "N's/100kbp":      f"{_st.n_per_100kbp:.1f}" if _st.n_per_100kbp is not None else "—",
                        "Core genes %":    _core_str,
                        "QC flags":        _flags,
                    })
                    _contig_paths[_sid] = Path(_cp)

        if _reads_derived_present:
            st.markdown("#### Assembly Results")
            _df_asm = pd.DataFrame(_asm_rows)
            if not _df_asm.empty:
                def _asm_row_style(row):
                    if str(row.get("Status", "")).startswith("❌"):
                        return ["background-color:#fef2f2;color:#dc2626"] * len(row)
                    if str(row.get("Status", "")).startswith("⚠️"):
                        return ["background-color:#fffbeb;color:#92400e"] * len(row)
                    return [""] * len(row)
                st.dataframe(
                    _df_asm.style.apply(_asm_row_style, axis=1),
                    use_container_width=True,
                    hide_index=True,
                )
            else:
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

