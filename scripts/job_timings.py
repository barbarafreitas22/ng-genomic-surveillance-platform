import sqlite3
import statistics
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from backend.paths import PROJECTS_DIR

JOBS_DB = PROJECTS_DIR / "jobs.db"


def _parse(ts: str | None) -> datetime | None:
    return datetime.fromisoformat(ts) if ts else None


def main() -> None:
    if not JOBS_DB.exists():
        print(f"No jobs database found at {JOBS_DB}")
        return

    con = sqlite3.connect(str(JOBS_DB))
    con.row_factory = sqlite3.Row
    rows = con.execute(
        "SELECT job_type, status, created_at, started_at, finished_at FROM jobs"
    ).fetchall()
    con.close()

    by_type: dict[str, list[dict]] = {}
    for r in rows:
        created  = _parse(r["created_at"])
        started  = _parse(r["started_at"])
        finished = _parse(r["finished_at"])
        entry = {
            "status":     r["status"],
            "queue_s":    (started - created).total_seconds() if created and started else None,
            "run_s":      (finished - started).total_seconds() if started and finished else None,
        }
        by_type.setdefault(r["job_type"], []).append(entry)

    if not by_type:
        print("No jobs recorded yet.")
        return

    for job_type, entries in sorted(by_type.items()):
        done   = [e for e in entries if e["status"] == "done"]
        errors = [e for e in entries if e["status"] == "error"]
        run_times = [e["run_s"] for e in done if e["run_s"] is not None]
        queue_times = [e["queue_s"] for e in done if e["queue_s"] is not None]

        print(f"\n=== {job_type} ===")
        print(f"  total: {len(entries)}  done: {len(done)}  error: {len(errors)}")

        if run_times:
            print(
                f"  run time (s)   — mean: {statistics.mean(run_times):.1f}  "
                f"median: {statistics.median(run_times):.1f}  "
                f"min: {min(run_times):.1f}  max: {max(run_times):.1f}"
            )
        else:
            print("  run time (s)   — no completed jobs")

        if queue_times:
            print(
                f"  queue wait (s) — mean: {statistics.mean(queue_times):.1f}  "
                f"median: {statistics.median(queue_times):.1f}  "
                f"min: {min(queue_times):.1f}  max: {max(queue_times):.1f}"
            )


if __name__ == "__main__":
    main()
