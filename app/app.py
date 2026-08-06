from __future__ import annotations

import sys
from pathlib import Path

import streamlit as st

APP_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = APP_DIR.parent
sys.path.insert(0, str(APP_DIR))
sys.path.append(str(PROJECT_ROOT))

try:
    from backend.logger import configure_logging as _configure_logging
    _configure_logging()
except Exception:
    pass

LOGO_PATH = str(APP_DIR / "logo.png")

st.set_page_config(
    page_title="Neisseria gonorrhoeae Genomic Surveillance Platform",
    page_icon=LOGO_PATH,
    layout="wide",
    initial_sidebar_state="expanded",
)

if not st.session_state.get("_entered_app"):
    st.markdown("""
    <style>
    @keyframes ngFadeUp {
        from { opacity: 0; transform: translateY(16px); }
        to   { opacity: 1; transform: translateY(0); }
    }
    [data-testid="stImage"] {
        display: flex;
        justify-content: center;
        width: 100%;
    }
    [data-testid="stImage"] img {
        animation: ngFadeUp 0.6s ease-out;
        transition: transform 0.25s ease;
    }
    [data-testid="stImage"] img:hover {
        transform: scale(1.06);
    }
    h2 {
        animation: ngFadeUp 0.6s ease-out 0.15s both;
    }
    </style>
    """, unsafe_allow_html=True)

    @st.cache_data
    def _reference_genome_stats() -> str | None:
        try:
            from Bio import SeqIO
            from backend.paths import REFERENCE_FA1090
            record = next(SeqIO.parse(str(REFERENCE_FA1090), "fasta"))
            seq = str(record.seq).upper()
            gc = (seq.count("G") + seq.count("C")) / len(seq) * 100
            return f"Reference genome: {record.id} · {len(seq):,} bp · {gc:.1f}% GC"
        except Exception:
            return None

    _c1, _c2, _c3 = st.columns([1, 2, 1])
    with _c2:
        st.image(LOGO_PATH, width=340)
        st.markdown(
            "<h2 style='text-align:center;'><em>Neisseria gonorrhoeae</em><br>Genomic Surveillance Platform</h2>",
            unsafe_allow_html=True,
        )
        _ref_stats = _reference_genome_stats()
        if _ref_stats:
            st.markdown(
                f"<p style='text-align:center;color:#64748b;font-size:0.85rem;'>{_ref_stats}</p>",
                unsafe_allow_html=True,
            )
        if st.button("Enter platform", type="primary", width="stretch"):
            st.session_state["_entered_app"] = True
            st.rerun()
        st.markdown(
            "<p style='text-align:center;margin-top:0.8rem;'>"
            "<a href='https://github.com/barbarafreitas22/ng-genomic-surveillance-platform' "
            "target='_blank' style='color:#64748b;text-decoration:none;font-size:0.85rem;'>"
            "&#128279; View on GitHub</a></p>",
            unsafe_allow_html=True,
        )
    st.stop()

st.sidebar.image(LOGO_PATH, width=140)

st.markdown("""
<style>
:root {
  --bg:          #ffffff;
  --bg-mid:      #f1f5f9;
  --bg-light:    #f8fafc;
  --border:      #e2e8f0;
  --teal:        #64748b;
  --teal-bright: #475569;
  --teal-dim:    #475569;
  --text-strong: #0f172a;
  --text-body:   #334155;
  --text-muted:  #64748b;
  --white:       #ffffff;
  --gray-50:     #f8fafc;
  --gray-100:    #f1f5f9;
  --gray-200:    #e2e8f0;
  --gray-300:    #cbd5e1;
  --gray-500:    #64748b;
  --gray-700:    #334155;
  --gray-900:    #0f172a;
  --radius-sm:   8px;
  --radius-md:   12px;
  --radius-lg:   16px;
  --shadow-sm:   0 1px 3px rgba(0,0,0,0.06), 0 1px 2px rgba(0,0,0,0.04);
  --shadow-md:   0 4px 12px rgba(0,0,0,0.08), 0 2px 4px rgba(0,0,0,0.04);
}

.block-container {
  padding-top: 1.5rem;
  padding-bottom: 3rem;
  max-width: 1200px;
}

body,
[data-testid="stAppViewContainer"],
[data-testid="stMain"],
[data-testid="stMain"] > div,
[data-testid="stHeader"],
[data-testid="stToolbar"] {
  background-color: var(--bg) !important;
}

p, span, li, td, th, label, div {
  color: var(--text-body);
}
h1, h2, h3, h4, h5, h6 { color: var(--text-strong) !important; }
h1 em, h2 em, h3 em { font-style: italic; }
.stMarkdown p, .stMarkdown li { color: var(--text-muted); }
code { color: var(--teal-dim) !important; background: rgba(100,116,139,0.08) !important; }

[data-testid="stTextInput"] input,
[data-testid="stSelectbox"] > div > div,
[data-testid="stTextArea"] textarea,
[data-testid="stNumberInput"] input {
  background: var(--white) !important;
  border-color: var(--border) !important;
  color: var(--text-strong) !important;
}
[data-testid="stTextInput"] input::placeholder { color: var(--text-muted) !important; }

[data-testid="stRadio"] label,
[data-testid="stCheckbox"] label { color: var(--text-body) !important; }

[data-testid="stCaptionContainer"] p { color: var(--text-muted) !important; }

hr { border-color: var(--border) !important; }

[data-testid="stAlert"] { border-radius: var(--radius-md) !important; }

[data-testid="stTooltipIcon"] { color: var(--text-muted) !important; }

[data-testid="stSidebar"] > div:first-child,
[data-testid="stSidebar"] {
  background: var(--bg) !important;
}
[data-testid="stSidebar"] .stMarkdown p,
[data-testid="stSidebar"] label,
[data-testid="stSidebar"] .stSelectbox label {
  color: var(--text-muted) !important;
  font-size: 0.73rem !important;
  font-weight: 600 !important;
  letter-spacing: 0.08em !important;
  text-transform: uppercase !important;
}
[data-testid="stSidebar"] h1 {
  color: var(--text-strong) !important;
  font-size: 0.8rem !important;
  font-weight: 700 !important;
  letter-spacing: 0.1em !important;
  text-transform: uppercase !important;
}
[data-testid="stSidebar"] .stSelectbox > div > div {
  background: var(--white) !important;
  border-color: var(--border) !important;
  color: var(--text-strong) !important;
}
[data-testid="stSidebar"] .stAlert {
  background: rgba(239,68,68,0.08) !important;
  border-color: rgba(239,68,68,0.25) !important;
}
[data-testid="stSidebar"] .stAlert p {
  color: #dc2626 !important;
  text-transform: none !important;
  font-size: 0.84rem !important;
  letter-spacing: 0 !important;
}
[data-testid="stSidebarNav"] { display: none; }

[data-testid="stFileUploadDropzone"],
[data-testid="stFileUploadDropzone"] > div {
  background: var(--bg-mid) !important;
  border-color: var(--border) !important;
}
[data-testid="stFileUploadDropzone"] small,
[data-testid="stFileUploadDropzone"] span,
[data-testid="stFileUploadDropzone"] p {
  color: var(--text-muted) !important;
}

.ag-root-wrapper,
.ag-header,
.ag-header-row,
.ag-header-cell,
.ag-row,
.ag-cell {
  background-color: var(--bg) !important;
  color: var(--text-body) !important;
}
.ag-header, .ag-header-row {
  border-bottom: 1px solid var(--border) !important;
}
.ag-row { border-bottom: 1px solid var(--border) !important; }
.ag-row-odd { background-color: var(--bg-mid) !important; }

button[data-testid="baseButton-primary"],
.stButton > button[kind="primary"] {
  background: var(--teal) !important;
  color: #ffffff !important;
  border: none !important;
  font-weight: 700 !important;
  border-radius: var(--radius-sm) !important;
  letter-spacing: 0.01em !important;
  box-shadow: 0 2px 8px rgba(100,116,139,0.22) !important;
  transition: all 0.15s ease !important;
}
button[data-testid="baseButton-primary"]:hover,
.stButton > button[kind="primary"]:hover {
  background: var(--teal-bright) !important;
  box-shadow: 0 4px 16px rgba(100,116,139,0.32) !important;
  transform: translateY(-1px) !important;
}

button[data-testid="baseButton-secondary"],
.stButton > button:not([kind="primary"]) {
  border: 1.5px solid var(--teal) !important;
  border-radius: var(--radius-sm) !important;
  background: rgba(100,116,139,0.07) !important;
  color: var(--teal-dim) !important;
  font-weight: 600 !important;
  transition: all 0.15s ease !important;
}
button[data-testid="baseButton-secondary"]:hover,
.stButton > button:not([kind="primary"]):hover {
  background: rgba(100,116,139,0.14) !important;
  border-color: var(--teal) !important;
  color: var(--teal-dim) !important;
  box-shadow: 0 2px 8px rgba(100,116,139,0.18) !important;
}

.stDownloadButton > button {
  border: 1.5px solid var(--teal) !important;
  color: var(--teal-dim) !important;
  border-radius: var(--radius-sm) !important;
  font-weight: 600 !important;
  background: rgba(100,116,139,0.07) !important;
  transition: all 0.15s ease !important;
}
.stDownloadButton > button:hover {
  background: rgba(100,116,139,0.14) !important;
  border-color: var(--teal) !important;
  color: var(--teal-dim) !important;
  box-shadow: 0 2px 8px rgba(100,116,139,0.18) !important;
}

[data-testid="stMetric"] {
  background: var(--bg-mid);
  border: 1px solid var(--border);
  border-radius: var(--radius-md);
  padding: 1.1rem 1.2rem;
  box-shadow: var(--shadow-sm);
  transition: box-shadow 0.15s;
}
[data-testid="stMetric"]:hover {
  box-shadow: var(--shadow-md);
  border-color: var(--teal) !important;
}
[data-testid="stMetricLabel"] > div {
  font-size: 0.7rem !important;
  font-weight: 700 !important;
  text-transform: uppercase !important;
  letter-spacing: 0.09em !important;
  color: var(--text-muted) !important;
}
[data-testid="stMetricValue"] > div {
  font-size: 1.85rem !important;
  font-weight: 800 !important;
  color: var(--text-strong) !important;
}

[data-testid="stExpander"] {
  background: var(--bg-mid) !important;
  border: 1px solid var(--border) !important;
  border-radius: var(--radius-md) !important;
  box-shadow: var(--shadow-sm) !important;
  overflow: hidden !important;
}
[data-testid="stExpander"] summary {
  font-weight: 600 !important;
  color: var(--text-body) !important;
}
[data-testid="stExpander"] summary:hover {
  color: var(--teal-dim) !important;
}

[data-testid="stDataFrame"] {
  border: 1px solid var(--border) !important;
  border-radius: var(--radius-md) !important;
  overflow: hidden !important;
  box-shadow: var(--shadow-sm) !important;
}

[data-testid="stFileUploadDropzone"] {
  border: 2px dashed var(--border) !important;
  border-radius: var(--radius-md) !important;
  background: var(--bg-mid) !important;
  transition: border-color 0.15s !important;
}
[data-testid="stFileUploadDropzone"]:hover {
  border-color: var(--teal) !important;
  background: rgba(100,116,139,0.04) !important;
}

[data-testid="stProgress"] > div > div > div {
  background: linear-gradient(90deg, var(--teal), var(--teal-bright)) !important;
}

.stTabs [data-baseweb="tab-list"] {
  border-bottom: 2px solid var(--border);
}
.stTabs [data-baseweb="tab"] {
  font-weight: 600 !important;
}
.stTabs [aria-selected="true"] {
  color: var(--teal-dim) !important;
}
.stTabs [data-baseweb="tab-highlight"] {
  background: var(--teal) !important;
}

hr {
  border-color: var(--border) !important;
}

.small-muted { color: var(--text-muted); font-size: 0.88rem; }
.ng-badge {
  display: inline-block;
  padding: 0.15rem 0.65rem;
  border-radius: 20px;
  font-size: 0.72rem;
  font-weight: 700;
  letter-spacing: 0.06em;
  text-transform: uppercase;
}
</style>
""", unsafe_allow_html=True)

from views.shared import (
    read_config, CONFIG_PATH, db,
    _load_project_cached, _rebuild_pipeline_results,
    ensure_result_page_state, init_directories,
)

ensure_result_page_state()
init_directories()

config = read_config(CONFIG_PATH)
paths = config["PATHS"] if config.has_section("PATHS") else {}
KRAKEN_DB_DEFAULT = paths.get("KRAKEN2_DB", "") if hasattr(paths, "get") else ""

try:
    _qp_project = st.query_params.get("project")
    _qp_page    = st.query_params.get("page")
    if _qp_project and db is not None:
        if st.session_state.get("active_project") != _qp_project:
            _qp_data = _load_project_cached(_qp_project)
            if _qp_data.get("samples"):
                st.session_state["active_project"]       = _qp_project
                st.session_state["full_pipeline_results"] = _rebuild_pipeline_results(_qp_data)
                if _qp_data.get("phylogeny"):
                    _qp_phy = _qp_data["phylogeny"]
                    st.session_state["last_tree_path"]      = _qp_phy.get("tree_path")
                    st.session_state["last_cluster_report"] = _qp_phy.get("cluster_report", {})
                    st.session_state["last_matrix_path"]    = _qp_phy.get("matrix_path")
                    _qp_cr = _qp_phy.get("cluster_report", {})
                    st.session_state["last_clusters"] = {
                        n: i["cluster"] for n, i in _qp_cr.items() if "cluster" in i
                    }
    if _qp_page == "results" and st.session_state.get("full_pipeline_results"):
        st.session_state["current_page"] = "Full Pipeline Results"
except Exception:
    pass

st.sidebar.markdown("""
<div style="padding:0.6rem 0 0.3rem 0;border-bottom:1px solid rgba(0,0,0,0.08);margin-bottom:0.8rem;">
  <div style="font-size:0.65rem;font-weight:700;letter-spacing:0.14em;text-transform:uppercase;
              color:#64748b;margin-bottom:0.15rem;">Navigation</div>
</div>
""", unsafe_allow_html=True)

pages_list = [
    "Homepage",
    "Full Pipeline Results",
    "1. Quality Control & Assembly",
    "2. AMR Profiling",
    "3. Phylogenetic Analysis",
    "4. Sample Metadata",
    "5. Neisseria gonorrhoeae Clinical Relevance",
    "6. Platform Technical Documentation",
]

if "_nav_to" in st.session_state:
    st.session_state["current_page"] = st.session_state.pop("_nav_to")

if st.session_state.get("current_page") not in pages_list:
    st.session_state["current_page"] = "Homepage"

page = st.sidebar.selectbox(
    "Select analysis module",
    pages_list,
    key="current_page",
)

_sidebar_proj = st.session_state.get("active_project")
if db and _sidebar_proj:
    try:
        _n_unread = db.count_unread_alerts(_sidebar_proj)
        if _n_unread > 0:
            _sev_icon = "🔴" if _n_unread else ""
            st.sidebar.error(f"{_sev_icon} **{_n_unread} unread alert(s)**")
    except Exception:
        pass

if page == "Homepage":
    from views.homepage import render; render()
elif page == "Full Pipeline Results":
    from views.full_pipeline_results import render; render()
elif page == "1. Quality Control & Assembly":
    from views.qc_assembly import render; render()
elif page == "2. AMR Profiling":
    from views.amr import render; render()
elif page == "3. Phylogenetic Analysis":
    from views.phylogenetic import render; render()
elif page == "4. Sample Metadata":
    from views.sample_metadata import render; render()
elif page == "5. Neisseria gonorrhoeae Clinical Relevance":
    from views.cr import render; render()
elif page == "6. Platform Technical Documentation":
    from views.documentation import render; render()
