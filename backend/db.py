import json
import sqlite3
import uuid
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path

from .paths import PROJECTS_DIR


def _db_path(project_name: str) -> Path:
    return PROJECTS_DIR / project_name / "project.db"


def _now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _migrate(con: sqlite3.Connection) -> None:
    new_cols = [
        ("qc_results",       "mean_coverage",           "REAL"),
        ("qc_results",       "pct_breadth_10x",         "REAL"),
        ("assembly_results", "misassemblies",            "INTEGER"),
        ("assembly_results", "duplication_ratio",        "REAL"),
        ("assembly_results", "nga50",                    "INTEGER"),
        ("assembly_results", "mismatches_per_100kbp",    "REAL"),
        ("assembly_results", "indels_per_100kbp",        "REAL"),
        ("assembly_results", "assembly_qc_status",       "TEXT"),
        ("assembly_results", "assembly_qc_flags",        "TEXT"),
        ("assembly_results", "plasmid_contigs_path",     "TEXT"),
        ("amr_results",      "ngstar_st",                "TEXT"),
        ("amr_results",      "ngstar_alleles",           "TEXT"),
        ("amr_results",      "ngstar_novel",             "INTEGER"),
        ("amr_results",      "ngstar_incomplete",        "INTEGER"),
        ("amr_results",      "resistance_category",      "TEXT"),
        ("amr_results",      "n_resistance_classes",     "INTEGER"),
        ("amr_results",      "essential_gene_mutations",  "TEXT"),
        ("amr_results",      "essential_gene_synonymous", "TEXT"),
        ("phylogeny_runs",   "cgmlst_results",           "TEXT"),
        ("phylogeny_runs",   "snp_clusters",             "TEXT"),
    ]
    for table, col, typ in new_cols:
        try:
            con.execute(f"ALTER TABLE {table} ADD COLUMN {col} {typ}")
        except sqlite3.OperationalError:
            pass


@contextmanager
def _conn(project_name: str):
    db = _db_path(project_name)
    db.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(str(db))
    con.row_factory = sqlite3.Row
    try:
        _migrate(con)
        yield con
        con.commit()
    finally:
        con.close()


def init_project(project_name: str) -> None:
    with _conn(project_name) as con:
        con.executescript("""
            CREATE TABLE IF NOT EXISTS samples (
                sample_id    TEXT PRIMARY KEY,
                added_at     TEXT NOT NULL,
                has_qc       INTEGER DEFAULT 0,
                has_assembly INTEGER DEFAULT 0,
                has_amr      INTEGER DEFAULT 0
            );

            CREATE TABLE IF NOT EXISTS qc_results (
                sample_id        TEXT PRIMARY KEY,
                total_reads      INTEGER,
                gc_content       INTEGER,
                trimmed_reads    INTEGER,
                ng_pct           REAL,
                neisseria_pct    REAL,
                unclassified_pct REAL,
                species_status   TEXT,
                mean_coverage    REAL,
                pct_breadth_10x  REAL,
                multiqc_path     TEXT,
                ran_at           TEXT
            );

            CREATE TABLE IF NOT EXISTS assembly_results (
                sample_id              TEXT PRIMARY KEY,
                contigs_path           TEXT,
                plasmid_contigs_path   TEXT,
                n_contigs              INTEGER,
                total_len              INTEGER,
                n50                    INTEGER,
                largest                INTEGER,
                gc_pct                 REAL,
                completeness           REAL,
                contigs_500            INTEGER,
                misassemblies          INTEGER,
                duplication_ratio      REAL,
                nga50                  INTEGER,
                mismatches_per_100kbp  REAL,
                indels_per_100kbp      REAL,
                assembly_qc_status     TEXT,
                assembly_qc_flags      TEXT,
                ran_at                 TEXT
            );

            CREATE TABLE IF NOT EXISTS amr_results (
                sample_id            TEXT PRIMARY KEY,
                failure_probability  REAL,
                cdc_phenotypes       TEXT,
                chromosomal          TEXT,
                plasmid              TEXT,
                therapy              TEXT,
                mlst_st              TEXT,
                mlst_alleles         TEXT,
                mosaic_pena_class    TEXT,
                mosaic_pena_identity REAL,
                mosaic_pena_coverage REAL,
                mosaic_suspected     INTEGER,
                ngstar_st            TEXT,
                ngstar_alleles       TEXT,
                ngstar_novel         INTEGER,
                ngstar_incomplete    INTEGER,
                resistance_category  TEXT,
                n_resistance_classes INTEGER,
                ran_at               TEXT
            );

            CREATE TABLE IF NOT EXISTS sample_metadata (
                sample_id       TEXT PRIMARY KEY,
                collection_date TEXT,
                anatomical_site TEXT,
                sex             TEXT,
                age             INTEGER,
                country         TEXT,
                city            TEXT,
                region          TEXT,
                health_unit     TEXT,
                updated_at      TEXT
            );

            CREATE TABLE IF NOT EXISTS phylogeny_runs (
                run_id         TEXT PRIMARY KEY,
                created_at     TEXT,
                tree_path      TEXT,
                matrix_path    TEXT,
                n_samples      INTEGER,
                cluster_report TEXT
            );

            CREATE TABLE IF NOT EXISTS phenotypic_amr (
                sample_id      TEXT NOT NULL,
                antibiotic     TEXT NOT NULL,
                mic_value      REAL,
                halo_mm        REAL,
                interpretation TEXT,
                updated_at     TEXT,
                PRIMARY KEY (sample_id, antibiotic)
            );

            CREATE TABLE IF NOT EXISTS alerts (
                id           INTEGER PRIMARY KEY AUTOINCREMENT,
                sample_id    TEXT NOT NULL,
                alert_type   TEXT NOT NULL,
                severity     TEXT NOT NULL,
                message      TEXT NOT NULL,
                created_at   TEXT NOT NULL,
                acknowledged INTEGER DEFAULT 0,
                UNIQUE(sample_id, alert_type)
            );
        """)


def save_qc(project_name: str, sample_id: str, qc_out: dict) -> None:
    met = qc_out.get("metrics", {})
    k2  = qc_out.get("kraken2", {})
    now = _now()
    with _conn(project_name) as con:
        con.execute("""
            INSERT INTO samples (sample_id, added_at, has_qc)
            VALUES (?, ?, 1)
            ON CONFLICT(sample_id) DO UPDATE SET has_qc=1
        """, (sample_id, now))
        con.execute("""
            INSERT OR REPLACE INTO qc_results
                (sample_id, total_reads, gc_content, trimmed_reads,
                 ng_pct, neisseria_pct, unclassified_pct, species_status,
                 mean_coverage, pct_breadth_10x, multiqc_path, ran_at)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
        """, (
            sample_id,
            met.get("total_reads"),
            met.get("gc"),
            met.get("trimmed_reads"),
            k2.get("ng_pct"),
            k2.get("neisseria_pct"),
            k2.get("unclassified_pct"),
            k2.get("species_status"),
            met.get("mean_coverage"),
            met.get("pct_breadth_10x"),
            qc_out.get("multiqc"),
            now,
        ))


def save_assembly(project_name: str, sample_id: str, contigs_path: str, stats) -> None:
    now = _now()
    _d  = stats.to_dict() if hasattr(stats, "to_dict") else (stats or {})
    qc  = _d.get("qc", {})
    with _conn(project_name) as con:
        con.execute("""
            INSERT INTO samples (sample_id, added_at, has_assembly)
            VALUES (?, ?, 1)
            ON CONFLICT(sample_id) DO UPDATE SET has_assembly=1
        """, (sample_id, now))
        con.execute("""
            INSERT OR REPLACE INTO assembly_results
                (sample_id, contigs_path,
                 n_contigs, total_len, n50,
                 largest, gc_pct, completeness, contigs_500,
                 misassemblies, duplication_ratio, nga50,
                 mismatches_per_100kbp, indels_per_100kbp,
                 assembly_qc_status, assembly_qc_flags, ran_at)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """, (
            sample_id, contigs_path,
            _d.get("n_contigs"), _d.get("total_len"), _d.get("n50"),
            _d.get("largest"), _d.get("gc_pct"),
            _d.get("completeness"), _d.get("contigs_500"),
            _d.get("misassemblies"), _d.get("duplication_ratio"),
            _d.get("nga50"), _d.get("mismatches_per_100kbp"),
            _d.get("indels_per_100kbp"),
            qc.get("status"), json.dumps(qc.get("flags", [])),
            now,
        ))


def save_amr(project_name: str, sample_id: str, amr_result) -> None:
    now    = _now()
    mlst   = amr_result.mlst
    mos    = amr_result.mosaic_pena
    ngstar = amr_result.ngstar
    with _conn(project_name) as con:
        con.execute("""
            INSERT INTO samples (sample_id, added_at, has_amr)
            VALUES (?, ?, 1)
            ON CONFLICT(sample_id) DO UPDATE SET has_amr=1
        """, (sample_id, now))
        con.execute("""
            INSERT OR REPLACE INTO amr_results
                (sample_id, failure_probability, cdc_phenotypes,
                 chromosomal, plasmid, therapy,
                 mlst_st, mlst_alleles,
                 mosaic_pena_class, mosaic_pena_identity,
                 mosaic_pena_coverage, mosaic_suspected,
                 ngstar_st, ngstar_alleles, ngstar_novel, ngstar_incomplete,
                 resistance_category, n_resistance_classes,
                 essential_gene_mutations, essential_gene_synonymous, ran_at)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """, (
            sample_id,
            amr_result.failure_probability,
            json.dumps(amr_result.cdc_phenotypes),
            json.dumps(amr_result.chromosomal),
            json.dumps(amr_result.plasmid),
            json.dumps(amr_result.therapy),
            mlst.get("st") if not mlst.get("error") else None,
            json.dumps(mlst.get("alleles", {})),
            mos.get("allele_class"),
            mos.get("identity"),
            mos.get("coverage_pct"),
            1 if mos.get("mosaic_suspected") else 0,
            ngstar.get("ST") if not ngstar.get("error") else None,
            json.dumps(ngstar.get("alleles", {})),
            1 if ngstar.get("novel") else 0,
            1 if ngstar.get("incomplete") else 0,
            amr_result.resistance_category,
            amr_result.n_resistance_classes,
            json.dumps(amr_result.essential_gene_mutations or {}),
            json.dumps(amr_result.essential_gene_synonymous or {}),
            now,
        ))


def delete_amr(project_name: str) -> None:
    with _conn(project_name) as con:
        con.execute("DELETE FROM amr_results")
        con.execute("UPDATE samples SET has_amr = 0")


def save_phylogeny(project_name: str, phy_out: dict) -> None:
    cgmlst = phy_out.get("cgmlst_result")
    snp    = phy_out.get("snp_clusters")
    with _conn(project_name) as con:
        con.execute("""
            INSERT OR REPLACE INTO phylogeny_runs
                (run_id, created_at, tree_path, matrix_path, n_samples,
                 cluster_report, cgmlst_results, snp_clusters)
            VALUES (?,?,?,?,?,?,?,?)
        """, (
            phy_out.get("run_id"),
            _now(),
            phy_out.get("tree_path"),
            phy_out.get("matrix_path"),
            phy_out.get("n_samples"),
            json.dumps(phy_out.get("cluster_report", {})),
            json.dumps(cgmlst) if cgmlst else None,
            json.dumps(snp)    if snp    else None,
        ))


def load_project(project_name: str) -> dict:
    db = _db_path(project_name)
    if not db.exists():
        return {"samples": {}, "phylogeny": None}

    with _conn(project_name) as con:
        samples = {}
        for row in con.execute("SELECT * FROM samples ORDER BY added_at"):
            sid = row["sample_id"]
            samples[sid] = {
                "added_at":     row["added_at"],
                "has_qc":       bool(row["has_qc"]),
                "has_assembly": bool(row["has_assembly"]),
                "has_amr":      bool(row["has_amr"]),
            }

        for row in con.execute("SELECT * FROM qc_results"):
            sid = row["sample_id"]
            if sid in samples:
                samples[sid]["qc"] = {
                    "metrics": {
                        "total_reads":     row["total_reads"],
                        "gc":              row["gc_content"],
                        "trimmed_reads":   row["trimmed_reads"],
                        "mean_coverage":   row["mean_coverage"],
                        "pct_breadth_10x": row["pct_breadth_10x"],
                    },
                    "kraken2": {
                        "ng_pct":           row["ng_pct"],
                        "neisseria_pct":    row["neisseria_pct"],
                        "unclassified_pct": row["unclassified_pct"],
                        "species_status":   row["species_status"],
                    },
                    "multiqc": row["multiqc_path"],
                }

        for row in con.execute("SELECT * FROM assembly_results"):
            sid = row["sample_id"]
            if sid in samples:
                samples[sid]["assembly"] = {
                    "contigs_path": row["contigs_path"],
                    "stats": {
                        "n_contigs":             row["n_contigs"],
                        "total_len":             row["total_len"],
                        "n50":                   row["n50"],
                        "largest":               row["largest"],
                        "gc_pct":                row["gc_pct"],
                        "completeness":          row["completeness"],
                        "contigs_500":           row["contigs_500"],
                        "misassemblies":         row["misassemblies"],
                        "duplication_ratio":     row["duplication_ratio"],
                        "nga50":                 row["nga50"],
                        "mismatches_per_100kbp": row["mismatches_per_100kbp"],
                        "indels_per_100kbp":     row["indels_per_100kbp"],
                        "qc": {
                            "status": row["assembly_qc_status"],
                            "flags":  json.loads(row["assembly_qc_flags"] or "[]"),
                        },
                    },
                }

        for row in con.execute("SELECT * FROM amr_results"):
            sid = row["sample_id"]
            if sid in samples:
                from .models import AMRResult as _AMRResult
                samples[sid]["amr"] = _AMRResult.from_dict({
                    "failure_probability":  row["failure_probability"],
                    "cdc_phenotypes":       json.loads(row["cdc_phenotypes"] or "[]"),
                    "chromosomal":          json.loads(row["chromosomal"] or "{}"),
                    "plasmid":              json.loads(row["plasmid"] or "{}"),
                    "therapy":              json.loads(row["therapy"] or "{}"),
                    "resistance_category":  row["resistance_category"],
                    "n_resistance_classes": row["n_resistance_classes"],
                    "mlst": {
                        "st":      row["mlst_st"],
                        "alleles": json.loads(row["mlst_alleles"] or "{}"),
                    },
                    "mosaic_pena": {
                        "allele_class":     row["mosaic_pena_class"],
                        "identity":         row["mosaic_pena_identity"],
                        "coverage_pct":     row["mosaic_pena_coverage"],
                        "mosaic_suspected": bool(row["mosaic_suspected"]),
                    },
                    "ngstar": {
                        "ST":         row["ngstar_st"],
                        "alleles":    json.loads(row["ngstar_alleles"] or "{}"),
                        "novel":      bool(row["ngstar_novel"]),
                        "incomplete": bool(row["ngstar_incomplete"]),
                    },
                    "essential_gene_mutations":  json.loads(row["essential_gene_mutations"] or "{}"),
                    "essential_gene_synonymous": json.loads(row["essential_gene_synonymous"] or "{}"),
                }, sample_id=sid)

        phylogeny = None
        row = con.execute(
            "SELECT * FROM phylogeny_runs ORDER BY created_at DESC LIMIT 1"
        ).fetchone()
        if row:
            _cgmlst_raw = row["cgmlst_results"] if "cgmlst_results" in row.keys() else None
            _snp_raw    = row["snp_clusters"]    if "snp_clusters"    in row.keys() else None
            phylogeny = {
                "run_id":         row["run_id"],
                "created_at":     row["created_at"],
                "tree_path":      row["tree_path"],
                "matrix_path":    row["matrix_path"],
                "n_samples":      row["n_samples"],
                "cluster_report": json.loads(row["cluster_report"] or "{}"),
                "cgmlst_result":  json.loads(_cgmlst_raw) if _cgmlst_raw else None,
                "snp_clusters":   json.loads(_snp_raw)    if _snp_raw    else None,
            }

    return {"samples": samples, "phylogeny": phylogeny}


_METADATA_FIELDS = (
    "collection_date", "anatomical_site", "sex", "age",
    "country", "city", "region", "health_unit",
)


def save_metadata(project_name: str, rows: list[dict]) -> None:
    now = _now()
    with _conn(project_name) as con:
        for row in rows:
            sid = row.get("sample_id")
            if not sid:
                continue
            con.execute("""
                INSERT INTO sample_metadata
                    (sample_id, collection_date, anatomical_site, sex, age,
                     country, city, region, health_unit, updated_at)
                VALUES (?,?,?,?,?,?,?,?,?,?)
                ON CONFLICT(sample_id) DO UPDATE SET
                    collection_date = excluded.collection_date,
                    anatomical_site = excluded.anatomical_site,
                    sex             = excluded.sex,
                    age             = excluded.age,
                    country         = excluded.country,
                    city            = excluded.city,
                    region          = excluded.region,
                    health_unit     = excluded.health_unit,
                    updated_at      = excluded.updated_at
            """, (
                sid,
                row.get("collection_date") or None,
                row.get("anatomical_site") or None,
                row.get("sex") or None,
                int(row["age"]) if str(row.get("age", "")).isdigit() else None,
                row.get("country") or None,
                row.get("city") or None,
                row.get("region") or None,
                row.get("health_unit") or None,
                now,
            ))


def load_metadata(project_name: str) -> dict[str, dict]:
    db = _db_path(project_name)
    if not db.exists():
        return {}
    with _conn(project_name) as con:
        try:
            return {
                row["sample_id"]: {f: row[f] for f in _METADATA_FIELDS}
                for row in con.execute("SELECT * FROM sample_metadata")
            }
        except Exception:
            return {}


def save_phenotypic_amr(project_name: str, sample_id: str, rows: list[dict]) -> None:
    now = _now()
    with _conn(project_name) as con:
        for row in rows:
            ab = row.get("antibiotic")
            if not ab:
                continue
            con.execute("""
                INSERT INTO phenotypic_amr
                    (sample_id, antibiotic, mic_value, halo_mm, interpretation, updated_at)
                VALUES (?,?,?,?,?,?)
                ON CONFLICT(sample_id, antibiotic) DO UPDATE SET
                    mic_value      = excluded.mic_value,
                    halo_mm        = excluded.halo_mm,
                    interpretation = excluded.interpretation,
                    updated_at     = excluded.updated_at
            """, (
                sample_id, ab,
                row.get("mic_value") if row.get("mic_value") not in (None, "") else None,
                row.get("halo_mm")   if row.get("halo_mm")   not in (None, "") else None,
                row.get("interpretation") or None,
                now,
            ))


def load_phenotypic_amr(project_name: str) -> dict[str, dict]:
    db = _db_path(project_name)
    if not db.exists():
        return {}
    with _conn(project_name) as con:
        try:
            result: dict = {}
            for row in con.execute("SELECT * FROM phenotypic_amr"):
                sid = row["sample_id"]
                if sid not in result:
                    result[sid] = {}
                result[sid][row["antibiotic"]] = {
                    "mic_value":      row["mic_value"],
                    "halo_mm":        row["halo_mm"],
                    "interpretation": row["interpretation"],
                }
            return result
        except Exception:
            return {}


def save_alerts(project_name: str, alerts: list[dict]) -> None:
    if not alerts:
        return
    now = _now()
    with _conn(project_name) as con:
        for a in alerts:
            con.execute("""
                INSERT INTO alerts
                    (sample_id, alert_type, severity, message, created_at)
                VALUES (?, ?, ?, ?, ?)
                ON CONFLICT(sample_id, alert_type) DO UPDATE SET
                    severity   = excluded.severity,
                    message    = excluded.message,
                    created_at = excluded.created_at
            """, (a["sample_id"], a["alert_type"], a["severity"], a["message"], now))


def load_alerts(project_name: str, unread_only: bool = False) -> list[dict]:
    db = _db_path(project_name)
    if not db.exists():
        return []
    with _conn(project_name) as con:
        try:
            where = "WHERE acknowledged=0" if unread_only else ""
            rows = con.execute(
                f"SELECT * FROM alerts {where} ORDER BY created_at DESC"
            ).fetchall()
            return [dict(r) for r in rows]
        except Exception:
            return []


def count_unread_alerts(project_name: str) -> int:
    db = _db_path(project_name)
    if not db.exists():
        return 0
    try:
        with _conn(project_name) as con:
            return con.execute(
                "SELECT COUNT(*) FROM alerts WHERE acknowledged=0"
            ).fetchone()[0]
    except Exception:
        return 0


def acknowledge_alerts(project_name: str, alert_ids: list[int] | None = None) -> None:
    with _conn(project_name) as con:
        if alert_ids is None:
            con.execute("UPDATE alerts SET acknowledged=1")
        else:
            con.executemany(
                "UPDATE alerts SET acknowledged=1 WHERE id=?",
                [(i,) for i in alert_ids],
            )


def list_projects() -> list[str]:
    if not PROJECTS_DIR.exists():
        return []
    return sorted(
        p.name for p in PROJECTS_DIR.iterdir()
        if p.is_dir() and (p / "project.db").exists()
    )


def projects_summary() -> list[dict]:
    rows = []
    for name in list_projects():
        try:
            with _conn(name) as con:
                n_qc  = con.execute("SELECT COUNT(*) FROM samples WHERE has_qc=1").fetchone()[0]
                n_asm = con.execute("SELECT COUNT(*) FROM samples WHERE has_assembly=1").fetchone()[0]
                n_amr = con.execute("SELECT COUNT(*) FROM samples WHERE has_amr=1").fetchone()[0]
                phy   = con.execute(
                    "SELECT n_samples, created_at FROM phylogeny_runs ORDER BY created_at DESC LIMIT 1"
                ).fetchone()
                rows.append({
                    "Project":   name,
                    "QC":        n_qc,
                    "Assembly":  n_asm,
                    "AMR":       n_amr,
                    "Phylogeny": f"{phy['n_samples']} samples" if phy else "—",
                    "Last tree": phy["created_at"][:10] if phy else "—",
                })
        except Exception:
            rows.append({
                "Project": name, "QC": 0, "Assembly": 0,
                "AMR": 0, "Phylogeny": "—", "Last tree": "—",
            })
    return rows


_JOBS_DB = PROJECTS_DIR / "jobs.db"


@contextmanager
def _jobs_conn():
    _JOBS_DB.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(str(_JOBS_DB), timeout=15)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("PRAGMA synchronous=NORMAL")
    try:
        yield con
        con.commit()
    finally:
        con.close()


def init_jobs_db(recover: bool = False) -> None:
    with _jobs_conn() as con:
        con.executescript("""
            CREATE TABLE IF NOT EXISTS jobs (
                id           TEXT PRIMARY KEY,
                project_name TEXT NOT NULL,
                sample_id    TEXT NOT NULL,
                job_type     TEXT NOT NULL,
                status       TEXT NOT NULL DEFAULT 'queued',
                payload      TEXT NOT NULL,
                result       TEXT,
                error_msg    TEXT,
                worker_id    TEXT,
                created_at   TEXT NOT NULL,
                started_at   TEXT,
                finished_at  TEXT
            );
            CREATE INDEX IF NOT EXISTS idx_jobs_status  ON jobs(status);
            CREATE INDEX IF NOT EXISTS idx_jobs_project ON jobs(project_name);
        """)
        if recover:
            con.execute(
                "UPDATE jobs SET status='queued', worker_id=NULL, started_at=NULL"
                " WHERE status='running'"
            )


def submit_job(project_name: str, sample_id: str, job_type: str, payload: dict) -> str:
    job_id = str(uuid.uuid4())
    init_jobs_db()
    with _jobs_conn() as con:
        con.execute(
            "INSERT INTO jobs (id,project_name,sample_id,job_type,payload,created_at) "
            "VALUES (?,?,?,?,?,?)",
            (job_id, project_name, sample_id, job_type, json.dumps(payload), _now()),
        )
    return job_id


def claim_next_job(worker_id: str) -> dict | None:
    now = _now()
    with _jobs_conn() as con:
        candidate = con.execute(
            "SELECT id FROM jobs WHERE status='queued' ORDER BY created_at LIMIT 1"
        ).fetchone()
        if candidate is None:
            return None
        job_id = candidate["id"]
        cur = con.execute(
            "UPDATE jobs SET status='running', worker_id=?, started_at=? "
            "WHERE id=? AND status='queued'",
            (worker_id, now, job_id),
        )
        if cur.rowcount == 0:
            return None
        row = con.execute("SELECT * FROM jobs WHERE id=?", (job_id,)).fetchone()
    return dict(row) if row else None


def complete_job(job_id: str, result: dict) -> None:
    with _jobs_conn() as con:
        con.execute(
            "UPDATE jobs SET status='done', result=?, finished_at=? WHERE id=?",
            (json.dumps(result), _now(), job_id),
        )


def fail_job(job_id: str, error: str) -> None:
    with _jobs_conn() as con:
        con.execute(
            "UPDATE jobs SET status='error', error_msg=?, finished_at=? WHERE id=?",
            (error, _now(), job_id),
        )


def get_project_jobs(project_name: str) -> list[dict]:
    try:
        init_jobs_db()
        with _jobs_conn() as con:
            rows = con.execute(
                "SELECT * FROM jobs WHERE project_name=? ORDER BY created_at",
                (project_name,),
            ).fetchall()
        return [dict(r) for r in rows]
    except Exception:
        return []


def cancel_job(job_id: str) -> None:
    with _jobs_conn() as con:
        con.execute(
            "UPDATE jobs SET status='error', error_msg='Cancelled by user' WHERE id=? AND status='queued'",
            (job_id,),
        )
