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
    _load_project_cached, _SITE_OPTIONS, _SEX_OPTIONS, _META_COLS,
)


def render() -> None:
    st.html("""
<div style="padding:0.6rem 0 1rem 0;">
  <div style="display:inline-flex;align-items:center;background:rgba(0,0,0,0.04);border:1px solid rgba(0,0,0,0.12);
  border-radius:20px;padding:0.22rem 0.85rem;font-size:0.7rem;font-weight:700;color:#64748b;
  letter-spacing:0.12em;text-transform:uppercase;margin-bottom:0.6rem;">Module 4</div>
  <div style="font-size:1.6rem;font-weight:800;color:#0f172a;">4. Sample Metadata</div>
</div>
""")

    _meta_proj = st.session_state.get("active_project")
    if not _meta_proj or not db:
        st.info("Open a project on the Home page first.")
    else:
        st.caption(f"Project: **{_meta_proj}**")

        _SITE_OPTIONS   = ["", "urethral", "rectal", "pharyngeal", "ocular", "conjunctival", "other"]
        _SEX_OPTIONS    = ["", "male", "female", "unknown"]
        _META_COLS      = ["sample_id", "collection_date", "anatomical_site", "sex", "age",
                           "country", "city", "region", "health_unit"]
        _META_LABELS    = {
            "sample_id": "Sample ID", "collection_date": "Collection date",
            "anatomical_site": "Anatomical site", "sex": "Sex", "age": "Age",
            "country": "Country", "city": "City/Municipality",
            "region": "Region", "health_unit": "Health unit",
        }

        db.init_project(_meta_proj)
        _proj_data   = _load_project_cached(_meta_proj)
        _known_sids  = sorted(_proj_data.get("samples", {}).keys())
        _existing    = db.load_metadata(_meta_proj)

        st.subheader("Import from CSV")
        st.caption(
            "The CSV must have a `sample_id` column. "
            "Optional columns: `collection_date` (YYYY-MM-DD), `anatomical_site`, "
            "`sex`, `age`, `country`, `city`, `region`, `health_unit`."
        )
        _meta_tmpl_rows = [",".join(_META_COLS)]
        if _known_sids:
            for _s in _known_sids[:3]:
                _meta_tmpl_rows.append(f"{_s}," + ",".join([""] * (len(_META_COLS) - 1)))
        st.download_button(
            "⬇ Download template CSV",
            data="\n".join(_meta_tmpl_rows).encode(),
            file_name="metadata_template.csv",
            mime="text/csv",
            key="meta_tmpl_dl",
        )

        _meta_csv = st.file_uploader(
            "Upload metadata CSV", type=["csv"], key="meta_csv_upload",
        )
        if _meta_csv is not None:
            try:
                _df_upload = pd.read_csv(_meta_csv, dtype=str).fillna("")
                if "sample_id" not in _df_upload.columns:
                    st.error("CSV must have a `sample_id` column.")
                else:
                    _valid_cols = [c for c in _META_COLS if c in _df_upload.columns]
                    _rows_to_save = _df_upload[_valid_cols].to_dict(orient="records")
                    db.save_metadata(_meta_proj, _rows_to_save)
                    st.success(f"{len(_rows_to_save)} sample(s) updated from CSV.")
                    _existing = db.load_metadata(_meta_proj)
            except Exception as _e:
                st.error(f"Error reading CSV: {_e}")

        st.divider()

        # ── Editable table ─────────────────────────────────────────────────
        st.subheader("Edit metadata")

        # Build dataframe: one row per known sample (fill blanks for missing)
        _blank = {f: "" for f in _META_COLS[1:]}
        _table_rows = []
        for _s in _known_sids:
            _row = {"sample_id": _s}
            _row.update({**_blank, **{k: (v or "") for k, v in _existing.get(_s, {}).items()}})
            _table_rows.append(_row)
        if not _table_rows:
            st.info("No samples found in this project. Run the pipeline first.")
        else:
            _df_edit = pd.DataFrame(_table_rows, columns=_META_COLS)
            _edited = st.data_editor(
                _df_edit,
                column_config={
                    "sample_id":       st.column_config.TextColumn("Sample ID", disabled=True),
                    "collection_date": st.column_config.TextColumn("Collection date", help="YYYY-MM-DD"),
                    "anatomical_site": st.column_config.SelectboxColumn("Anatomical site", options=_SITE_OPTIONS),
                    "sex":             st.column_config.SelectboxColumn("Sex", options=_SEX_OPTIONS),
                    "age":             st.column_config.NumberColumn("Age", min_value=0, max_value=120, step=1),
                    "country":         st.column_config.TextColumn("Country"),
                    "city":            st.column_config.TextColumn("City/Municipality"),
                    "region":          st.column_config.TextColumn("Region"),
                    "health_unit":     st.column_config.TextColumn("Health unit"),
                },
                use_container_width=True,
                hide_index=True,
                num_rows="fixed",
                key="meta_editor",
            )

            if st.button("Save changes", type="primary", key="meta_save_btn"):
                _save_rows = _edited.to_dict(orient="records")
                db.save_metadata(_meta_proj, _save_rows)
                st.success(f"Metadata saved for {len(_save_rows)} sample(s).")
                _existing = db.load_metadata(_meta_proj)

        if _existing:
            st.divider()
            st.subheader("Current metadata")
            _summary_rows = []
            for _s in _known_sids:
                _m = _existing.get(_s, {})
                _summary_rows.append({
                    "Sample":           _s,
                    "Date":             _m.get("collection_date") or "—",
                    "Site":             _m.get("anatomical_site") or "—",
                    "Sex":              _m.get("sex") or "—",
                    "Age":              str(_m.get("age") or "—"),
                    "Country":          _m.get("country") or "—",
                    "City":             _m.get("city") or "—",
                    "Region":           _m.get("region") or "—",
                    "Health unit":      _m.get("health_unit") or "—",
                })
            _df_summary = pd.DataFrame(_summary_rows)
            st.dataframe(_df_summary, use_container_width=True, hide_index=True)
            st.download_button(
                "⬇ Export metadata CSV",
                data=_df_summary.to_csv(index=False).encode(),
                file_name=f"{_meta_proj}_metadata.csv",
                mime="text/csv",
                key="meta_export_dl",
            )


