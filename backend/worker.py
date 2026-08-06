import json
import logging
import os
import signal
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

try:
    import psutil as _psutil
except ImportError:
    _psutil = None

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [worker] %(levelname)s %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
)
log = logging.getLogger(__name__)

WORKER_ID     = f"worker-{os.getpid()}-{uuid.uuid4().hex[:6]}"
POLL_INTERVAL = float(os.getenv("WORKER_POLL_INTERVAL", "2"))
MAX_THREADS   = int(os.getenv("WORKER_THREADS", max(1, (os.cpu_count() or 4) // 4)))
MIN_FREE_GB   = float(os.getenv("WORKER_MIN_FREE_GB", "1.5"))


def _available_memory_gb() -> float | None:
    try:
        _cg_limit   = int(open("/sys/fs/cgroup/memory.max").read().strip())
        _cg_current = int(open("/sys/fs/cgroup/memory.current").read().strip())
        return (_cg_limit - _cg_current) / 1e9
    except Exception:
        return (_psutil.virtual_memory().available / 1e9) if _psutil is not None else None


def _run_job(job: dict) -> dict:
    from backend.db import init_project, save_qc, save_assembly, save_amr

    payload = json.loads(job["payload"])
    jtype   = job["job_type"]
    project = job["project_name"]
    sid     = job["sample_id"]
    init_project(project)

    if jtype == "qc":
        from backend.qc import run_qc_pipeline
        result = run_qc_pipeline(
            sid,
            [Path(p) for p in payload["files"]],
            output_dir=payload.get("output_dir"),
        )
        save_qc(project, sid, result)
        return result

    if jtype == "assembly":
        from backend.assembly import run_assembly_pipeline, parse_contigs_stats
        asm = run_assembly_pipeline(
            sample_id=sid,
            r1=Path(payload["r1"]),
            r2=Path(payload["r2"]) if payload.get("r2") else None,
            run_plasmid=payload.get("run_plasmid", True),
            pre_trimmed=payload.get("pre_trimmed", False),
        )
        result = {
            "contigs_path":         str(asm.contigs) if asm.contigs else None,
            "plasmid_contigs_path": str(asm.plasmid_contigs) if asm.plasmid_contigs else None,
            "logs": asm.logs,
        }
        if asm.contigs:
            from backend.phylogeny.cgmlst import core_genome_completeness
            try:
                core_genome_completeness(str(asm.contigs), str(Path(asm.contigs).parent))
            except Exception:
                pass
            stats = parse_contigs_stats(asm.contigs)
            save_assembly(project, sid, str(asm.contigs), stats,
                          plasmid_contigs_path=str(asm.plasmid_contigs) if asm.plasmid_contigs else None)
        return result

    if jtype == "amr":
        from backend.amr import run_amr_variant_calling
        from backend.models import AMRResult
        result = run_amr_variant_calling(
            sample_id=sid,
            contigs=Path(payload["contigs"]),
            plasmid_contigs=Path(payload["plasmid_contigs"]) if payload.get("plasmid_contigs") else None,
        )
        if isinstance(result, AMRResult):
            save_amr(project, sid, result)
            return result.to_dict()
        return result

    if jtype == "qc_fasta":
        from backend.qc import run_qc_pipeline_fasta
        from backend.assembly import parse_contigs_stats
        from backend.amr import run_amr_variant_calling
        from backend.models import AMRResult
        fasta_path = Path(payload["fasta_path"])

        qc_result = run_qc_pipeline_fasta(sid, fasta_path, output_dir=payload.get("output_dir"))
        save_qc(project, sid, qc_result)

        from backend.phylogeny.cgmlst import core_genome_completeness
        try:
            core_genome_completeness(str(fasta_path), str(fasta_path.parent))
        except Exception:
            pass
        stats = parse_contigs_stats(fasta_path)
        stats_dict = stats.to_dict() if stats else {}
        save_assembly(project, sid, str(fasta_path), stats_dict)

        try:
            amr_out = run_amr_variant_calling(sample_id=sid, contigs=fasta_path)
            save_amr(project, sid, amr_out)
            amr_dict = amr_out.to_dict() if isinstance(amr_out, AMRResult) else amr_out
        except Exception as e:
            amr_dict = {"error": str(e)}

        return {"qc": qc_result, "assembly_stats": stats_dict, "amr": amr_dict}

    if jtype == "full_pipeline":
        from backend.full_pipeline import run_full_pipeline
        from backend.models import AMRResult
        result = run_full_pipeline(
            project_name=project,
            input_paths=[Path(p) for p in payload["input_paths"]],
            run_plasmid=payload.get("run_plasmid", False),
        )
        result["amr"] = {
            sid: (amr.to_dict() if isinstance(amr, AMRResult) else amr)
            for sid, amr in result["amr"].items()
        }
        return result

    if jtype == "phylogeny_nj":
        from backend.phylogeny.run_phylogeny import build_distance_tree
        return build_distance_tree(payload["samples"])

    raise ValueError(f"Unknown job type: {jtype!r}")


def _run_and_finish(job: dict) -> None:
    from backend.db import complete_job, fail_job
    log.info("Starting job %s  type=%-12s  sample=%s  project=%s",
             job["id"][:8], job["job_type"], job["sample_id"], job["project_name"])
    try:
        result = _run_job(job)
        complete_job(job["id"], result)
        log.info("Job %s done", job["id"][:8])
    except Exception as exc:
        log.exception("Job %s failed: %s", job["id"][:8], exc)
        fail_job(job["id"], str(exc))


def main() -> None:
    from backend.db import claim_next_job, init_jobs_db

    init_jobs_db(recover=True)
    log.info("Worker %s started — %d threads — polling every %.0fs",
             WORKER_ID, MAX_THREADS, POLL_INTERVAL)

    running = True

    def _stop(sig, _frame):
        nonlocal running
        running = False

    signal.signal(signal.SIGTERM, _stop)
    signal.signal(signal.SIGINT, _stop)

    with ThreadPoolExecutor(max_workers=MAX_THREADS) as pool:
        active: set = set()
        while running or active:
            active -= {f for f in active if f.done()}
            while running and len(active) < MAX_THREADS:
                _avail = _available_memory_gb()
                if _avail is not None and _avail < MIN_FREE_GB and active:
                    log.warning("Low memory (%.1f GB free < %.1f GB threshold) — "
                                "holding off new jobs until one finishes", _avail, MIN_FREE_GB)
                    break
                job = claim_next_job(WORKER_ID)
                if job is None:
                    break
                active.add(pool.submit(_run_and_finish, job))
            time.sleep(POLL_INTERVAL)

    log.info("Worker %s stopped", WORKER_ID)


if __name__ == "__main__":
    main()
