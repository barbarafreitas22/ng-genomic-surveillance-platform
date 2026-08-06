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
    _load_project_cached,
    _render_cohort_amr_profile, _inline_metadata_widget,
    _render_essential_gene_markers,
)
from backend.models import AMRResult as _AMRResult


def render() -> None:
    st.html("""
<div style="padding:0.6rem 0 1rem 0;">
  <div style="display:inline-flex;align-items:center;background:rgba(0,0,0,0.04);border:1px solid rgba(0,0,0,0.12);
  border-radius:20px;padding:0.22rem 0.85rem;font-size:0.7rem;font-weight:700;color:#64748b;
  letter-spacing:0.12em;text-transform:uppercase;margin-bottom:0.6rem;">Module 2</div>
  <div style="font-size:1.6rem;font-weight:800;color:#0f172a;">2. AMR Profiling</div>
</div>
""")

    _amr_proj = st.session_state.get("active_project")
    if db and _amr_proj and st.session_state.get("_amr_loaded_project") != _amr_proj:
        _amr_db = _load_project_cached(_amr_proj)
        _amr_from_db = {sid: d["amr"] for sid, d in _amr_db["samples"].items() if "amr" in d}
        st.session_state["amr_manual_results"] = _amr_from_db
        st.session_state["amr_manual_bytes"]   = {}
        st.session_state["amr_manual_errors"]  = {}
        st.session_state["_amr_loaded_project"] = _amr_proj

    if _amr_proj:
        st.caption(f"Project: **{_amr_proj}** — results saved automatically")
    else:
        st.caption("No project selected — results will not be saved. Open a project on the Home page first.")

    st.markdown("""
    This module screens assembled *N. gonorrhoeae* contigs for antimicrobial resistance determinants,
    combining two complementary sources of evidence:

    - **Chromosomal mutations** — point mutations in resistance-associated genes, reported in amino-acid
      notation (e.g. *gyrA* S91F)
    - **Plasmid-borne genes** — acquired resistance genes, reported as presence/absence (e.g. *blaTEM-1*)

    The pipeline aligns each assembly against a curated resistance-gene panel with **Minimap2**, calls
    variants with **BCFtools**, and interprets the results against the **European 2020 (IUSTI) treatment
    guidelines** to assign a resistance category and therapy recommendation.
    """)

    def render_amr_summary_table(results_dict: dict, bytes_map: dict | None = None,
                                  key_prefix: str = "amr", project_name: str | None = None):
        valid = {sid: r for sid, r in results_dict.items() if isinstance(r, _AMRResult)}
        errors = {sid: r["error"] for sid, r in results_dict.items() if isinstance(r, dict) and "error" in r}

        _sel_key = f"{key_prefix}_selected_sid"
        _active_sid = st.session_state.get(_sel_key)

        if _active_sid and _active_sid in valid:
            if st.button("← Back to overview", key=f"{key_prefix}_back_btn"):
                st.session_state[_sel_key] = None
                st.rerun()
            st.markdown(f"### {_active_sid}")
            _tbl_ctg = st.session_state.get("assembly_results", {}).get(_active_sid, {}).get("contigs_path")
            if bytes_map and _active_sid in bytes_map:
                fname, fbytes = bytes_map[_active_sid]
                _dc, _vc, _ = st.columns([2, 2, 5])
                _dc.download_button(
                    "⬇ Download genome", data=fbytes, file_name=fname,
                    mime="application/octet-stream", key=f"{key_prefix}_dl_{_active_sid}",
                )
                if _vc.button("View FASTA", key=f"{key_prefix}_view_{_active_sid}"):
                    _fk = f"{key_prefix}_show_fasta_{_active_sid}"
                    st.session_state[_fk] = not st.session_state.get(_fk, False)
                if st.session_state.get(f"{key_prefix}_show_fasta_{_active_sid}"):
                    lines = fbytes.decode("utf-8", errors="replace").splitlines()
                    preview = "\n".join(lines[:50])
                    if len(lines) > 50:
                        preview += f"\n… ({len(lines) - 50} more lines)"
                    st.code(preview, language=None)
                st.markdown("")
            render_amr_interpretation(valid[_active_sid], sample_label=_active_sid, contigs_path=_tbl_ctg)
            return

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
            "efflux_pump_overexpression":          "MtrCDE efflux — reduced susceptibility (penicillin, tetracycline, azithromycin)",
            "norM_efflux_upregulation":            "NorM efflux — reduced fluoroquinolone susceptibility",
            "efflux_tetracycline_contribution":    None,
            "penA_allele_likely_mosaic":           None,
            "fluoroquinolone_minor_parE":          None,
            "fluoroquinolone_minor_gyrB":          None,
            "zoliflodacin_resistance_gyrB":        None,
        }

        def _pheno_label(p: str) -> str | None:
            if p == "wildtype":
                return None
            if p in _PHENOTYPE_LABEL:
                return _PHENOTYPE_LABEL[p]
            return p.replace("_", " ").capitalize()

        summary_rows: list = []
        sample_order: list = []

        for sid, res in valid.items():
            prob        = res.failure_probability
            cdc         = res.cdc_phenotypes
            ther        = res.therapy
            mlst        = res.mlst
            ngstar      = res.ngstar
            who_matches = res.who_matches
            res_cat     = res.resistance_category

            labels     = [_pheno_label(p) for p in cdc]
            labels     = [l for l in labels if l is not None]
            phenotypes = "\n".join(labels) or "Wildtype"

            treatment  = (ther.get("recommend") or ["—"])[0] if ther else "—"
            st_val     = mlst.get("st") or "—" if not mlst.get("error") else "—"
            ngstar_val = f"ST-{ngstar['ST']}" if not ngstar.get("error") and ngstar.get("ST") else "—"

            has_xdr = any(m["mdr_class"] == "XDR" for m in who_matches)
            is_mdr_xdr = has_xdr or res_cat == "MDR"
            sample_label = ("🔴 " + sid) if is_mdr_xdr else sid

            summary_rows.append({
                "Sample":                     sample_label,
                "Genomic AMR Interpretation": treatment,
                "CDC phenotype":              phenotypes,
                "Score":                      f"{prob*100:.1f}%",
                "ST (MLST)":                  st_val,
                "NG-STAR":                    ngstar_val,
                "_prob":                      prob,
            })
            sample_order.append(sid)

        for sid in errors:
            summary_rows.append({
                "Sample":                     sid,
                "Genomic AMR Interpretation": errors[sid][:80],
                "CDC phenotype":              "—",
                "Score":                      "—",
                "ST (MLST)":                  "—",
                "NG-STAR":                    "—",
                "_prob":                      -1.0,
            })
            sample_order.append(sid)

        if not summary_rows:
            return

        render_amr_mutation_matrix(results_dict)
        _render_essential_gene_markers(results_dict)

        df_sum = pd.DataFrame(summary_rows)
        display_df = df_sum.drop(columns=["_prob"])

        def _tbl_row_style(row):
            if df_sum.loc[row.name, "_prob"] < 0:
                return ["background-color:#fef2f2;color:#dc2626"] * len(row)
            return [""] * len(row)

        event = st.dataframe(
            display_df.style.apply(_tbl_row_style, axis=1),
            use_container_width=True,
            hide_index=True,
            selection_mode="single-row",
            on_select="rerun",
            key=f"{key_prefix}_tbl",
        )
        st.download_button(
            "⬇ Download summary report (CSV)",
            data=display_df.to_csv(index=False).encode(),
            file_name="amr_summary.csv",
            mime="text/csv",
            key=f"{key_prefix}_csv_dl",
        )

        selected_rows = (event.selection.rows
                         if event and hasattr(event, "selection") else [])

        if selected_rows:
            sel_sid = sample_order[selected_rows[0]]
            if sel_sid in valid and st.session_state.get(_sel_key) != sel_sid:
                st.session_state[_sel_key] = sel_sid
                st.rerun()
        else:
            st.caption("Click a row in the table above to view the full report for that sample.")

    pipeline_results = st.session_state.get("full_pipeline_results")
    amr_pipeline = pipeline_results.get("amr", {}) if pipeline_results else {}

    if amr_pipeline:
        valid_pipe = {sid: r for sid, r in amr_pipeline.items() if isinstance(r, _AMRResult)}
        if len(valid_pipe) == 1:
            sid = next(iter(valid_pipe))
            _pipe_ctg = (pipeline_results or {}).get("assembly", {}).get(sid, {}).get("contigs_path")
            render_amr_interpretation(valid_pipe[sid], sample_label=sid, contigs_path=_pipe_ctg)
        elif valid_pipe:
            render_amr_summary_table(amr_pipeline, key_prefix="pipe", project_name=_amr_proj)
        else:
            for sid, r in amr_pipeline.items():
                st.error(f"**{sid}**: {r.get('error', 'unknown error')}")
        st.divider()

    #  Manual results
    manual_results = {
        sid: r for sid, r in st.session_state.get("amr_manual_results", {}).items()
        if sid not in amr_pipeline
    }
    manual_bytes   = st.session_state.get("amr_manual_bytes",   {})
    manual_errors  = st.session_state.get("amr_manual_errors",  {})

    if manual_errors:
        for sid_err, msg in manual_errors.items():
            st.error(f"**{sid_err}**: {msg}")

    if manual_results:
        if len(manual_results) == 1:
            sid = next(iter(manual_results))
            if sid in manual_bytes:
                fname, fbytes = manual_bytes[sid]
                _dc, _vc, _ = st.columns([2, 2, 5])
                _dc.download_button(
                    "⬇ Download genome",
                    data=fbytes,
                    file_name=fname,
                    mime="application/octet-stream",
                    key="man_dl_single",
                )
                if _vc.button("View FASTA", key="man_view_single"):
                    st.session_state["man_show_fasta_single"] = not st.session_state.get("man_show_fasta_single", False)
                if st.session_state.get("man_show_fasta_single"):
                    lines = fbytes.decode("utf-8", errors="replace").splitlines()
                    preview = "\n".join(lines[:50])
                    if len(lines) > 50:
                        preview += f"\n… ({len(lines) - 50} more lines)"
                    st.code(preview, language=None)
                st.markdown("")
            render_amr_interpretation(manual_results[sid], sample_label=sid)
        else:
            render_amr_summary_table(manual_results, bytes_map=manual_bytes, key_prefix="man", project_name=_amr_proj)
        st.divider()

    if _amr_proj and db:
        _amr_cohort_data = {
            sid: d["amr"]
            for sid, d in _load_project_cached(_amr_proj)["samples"].items()
            if "amr" in d
        }
        _render_cohort_amr_profile(_amr_cohort_data, project_name=_amr_proj)

    has_any_results = bool(amr_pipeline) or bool(manual_results)

    if has_any_results and _amr_proj and run_amr_variant_calling:
        st.divider()
        _rcol, _ = st.columns([3, 4])
        if _rcol.button("Delete results and re-run AMR for this project", key="amr_rerun_btn", type="primary"):
            # Get contigs paths from DB assembly results
            _proj_data = _load_project_cached(_amr_proj)
            _contigs_map = {
                sid: Path(d["assembly"]["contigs_path"])
                for sid, d in _proj_data["samples"].items()
                if "assembly" in d and d["assembly"].get("contigs_path")
            }
            _plasmid_map = {
                sid: Path(d["assembly"]["plasmid_contigs_path"])
                for sid, d in _proj_data["samples"].items()
                if "assembly" in d and d["assembly"].get("plasmid_contigs_path")
            }

            from backend.db import delete_amr as _delete_amr
            try:
                _delete_amr(_amr_proj)
            except Exception as _del_err:
                st.error(f"Failed to delete AMR results: {_del_err}")
                st.stop()

            if pipeline_results and "amr" in pipeline_results:
                pipeline_results["amr"] = {}
            st.session_state["amr_manual_results"] = {}
            st.session_state["amr_manual_bytes"]   = {}
            st.session_state["amr_manual_errors"]  = {}
            st.session_state["_amr_loaded_project"] = None
            _load_project_cached.clear()

            if not _contigs_map:
                st.info("Results deleted. No assembled contigs found — upload genomes via 'Run on new samples'.")
                st.rerun()

            _rerun_results: dict = {}
            _rerun_errors:  dict = {}
            _n_rr = len(_contigs_map)
            _workers_rr = max(1, min(_n_rr, (os.cpu_count() or 2) // 2))
            _prog = st.progress(0, text=f"Re-running AMR for {_n_rr} sample(s) — {_workers_rr} in parallel…")
            _done_rr = 0
            with ThreadPoolExecutor(max_workers=_workers_rr) as _exec_rr:
                _rr_futures = {
                    _exec_rr.submit(
                        run_amr_variant_calling, _sid, _cpath,
                        _plasmid_map.get(_sid),
                    ): _sid
                    for _sid, _cpath in _contigs_map.items()
                }
                for _fut in as_completed(_rr_futures):
                    _sid = _rr_futures[_fut]
                    _done_rr += 1
                    _prog.progress(_done_rr / _n_rr, text=f"Done {_done_rr}/{_n_rr} — {_sid}")
                    try:
                        _res = _fut.result()
                        _rerun_results[_sid] = _res
                    except Exception as _exc:
                        _rerun_errors[_sid] = str(_exc)
            _prog.progress(1.0, text="Done")

            if db:
                for _sid, _res in _rerun_results.items():
                    if isinstance(_res, _AMRResult):
                        db.save_amr(_amr_proj, _sid, _res)
                _load_project_cached.clear()
            st.session_state["amr_manual_results"] = _rerun_results
            st.session_state["amr_manual_errors"]  = _rerun_errors
            st.session_state["_amr_loaded_project"] = None
            _n_ok  = len(_rerun_results)
            _n_err = len(_rerun_errors)
            if _n_err:
                st.warning(f"Re-run complete: {_n_ok} succeeded, {_n_err} failed.")
            else:
                st.success(f"Re-run complete — {_n_ok} sample(s) analysed.")
            st.rerun()

    # Upload
    with st.expander("Run on new samples", expanded=not has_any_results):
        contigs_files = st.file_uploader(
            "Upload one or more assembled genomes (contigs.fasta)",
            type=["fasta", "fa", "fna"],
            accept_multiple_files=True,
            key="amr_multi_upload",
        )

        amr_id_map: dict = {}
        amr_bytes_map: dict = {}

        if contigs_files:
            st.markdown("**Sample IDs** — edit if needed before running:")
            for i, f in enumerate(contigs_files):
                c1, c2 = st.columns([3, 2])
                c1.markdown(f"`{f.name}`  ({f.size / 1024:.1f} KB)")
                default_id = Path(f.name).stem
                sid_input = c2.text_input(
                    "ID",
                    value=default_id,
                    placeholder="e.g. NG-001",
                    key=f"amr_sid_{i}",
                    label_visibility="collapsed",
                )
                sid_clean = sanitize_name(sid_input) if sid_input.strip() else default_id
                amr_id_map[f.name] = sid_clean
                amr_bytes_map[sid_clean] = (f.name, bytes(f.getbuffer()))

            _inline_metadata_widget(
                list(amr_id_map.values()), _amr_proj, "amr_up"
            )

            if st.button("Run AMR Analysis", type="primary", key="amr_run_btn"):
                run_results: dict = {}
                run_errors: dict = {}
                _n_amr = len(contigs_files)
                _workers_amr = max(1, min(_n_amr, (os.cpu_count() or 2) // 2))
                prog = st.progress(0, text="Saving uploads…")

                # Save all uploads on the main thread before parallelising
                _amr_saved_map: dict[str, Path] = {}
                _amr_all_saved: list[Path] = []
                for f in contigs_files:
                    sid = amr_id_map[f.name]
                    _path = save_temp_file(f)
                    _amr_saved_map[sid] = _path
                    _amr_all_saved.append(_path)

                prog.progress(0, text=f"Analysing {_n_amr} sample(s) — {_workers_amr} in parallel…")
                _done_amr = 0
                try:
                    with ThreadPoolExecutor(max_workers=_workers_amr) as _exec_amr:
                        _amr_futures = {
                            _exec_amr.submit(run_amr_variant_calling, _sid, _path): _sid
                            for _sid, _path in _amr_saved_map.items()
                        }
                        for _fut in as_completed(_amr_futures):
                            _sid = _amr_futures[_fut]
                            _done_amr += 1
                            prog.progress(_done_amr / _n_amr, text=f"Done {_done_amr}/{_n_amr} — {_sid}")
                            try:
                                _res = _fut.result()
                                run_results[_sid] = _res
                            except Exception as exc:
                                run_errors[_sid] = str(exc)
                    prog.progress(1.0, text="Done")
                finally:
                    cleanup_temp_files(_amr_all_saved)

                st.session_state["amr_manual_results"] = run_results
                st.session_state["amr_manual_bytes"]   = amr_bytes_map
                st.session_state["amr_manual_errors"]  = run_errors
                if db and _amr_proj:
                    db.init_project(_amr_proj)
                    for sid, res in run_results.items():
                        if isinstance(res, _AMRResult):
                            db.save_amr(_amr_proj, sid, res)
                    st.session_state["_amr_loaded_project"] = _amr_proj
                    _load_project_cached.clear()
                st.rerun()


