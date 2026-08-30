import logging
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import re

logger = logging.getLogger(__name__)

try:
    import psutil as _psutil
except ImportError:
    _psutil = None

from backend.db import init_project, save_qc, save_assembly, save_amr, save_phylogeny
from backend.qc import run_qc_pipeline, run_qc_pipeline_fasta
from backend.assembly import run_assembly_pipeline, parse_contigs_stats as _parse_contigs_stats
from backend.amr import run_amr_variant_calling
from backend.phylogeny.run_phylogeny import build_distance_tree
from backend.models import AMRResult

from .paths import PROJECTS_DIR as _PROJECTS_DIR


def _project_results_dir(project_name: str) -> Path:
    return _PROJECTS_DIR / project_name / "results"

def project_qc_dir(project_name: str) -> Path:
    return _project_results_dir(project_name) / "qc"

def project_assembly_dir(project_name: str) -> Path:
    return _project_results_dir(project_name) / "assembly"

def project_amr_dir(project_name: str) -> Path:
    return _project_results_dir(project_name) / "amr"

def project_phylogeny_dir(project_name: str) -> Path:
    return _project_results_dir(project_name) / "phylogeny"


def _detect_read_type(filename: str) -> str:
    name = filename.lower()
    if re.search(r"(?:^|[_.-])(r?1)(?:[_.-]|\.|$)", name):
        return "R1"
    if re.search(r"(?:^|[_.-])(r?2)(?:[_.-]|\.|$)", name):
        return "R2"
    return "SE"




def _group_samples(input_paths: list[Path]) -> dict[str, dict]:
    """Group files into samples, detecting R1/R2 pairs automatically."""
    samples: dict[str, dict] = {}
    for path in input_paths:
        rtype = _detect_read_type(path.name)
        base = re.sub(r"[_.-]?(?:r?[12])(?:\.[^.]+)*$", "", path.stem, flags=re.IGNORECASE)
        base = base or path.stem
        if base not in samples:
            samples[base] = {"R1": None, "R2": None}
        if rtype == "R2":
            samples[base]["R2"] = path
        else:
            samples[base]["R1"] = path
    return samples


def _process_single_sample(project_name: str, sample_id: str, r1: Path, r2) -> dict:
    result: dict = {"qc": None, "assembly": None, "amr": None, "contigs_path": None, "excluded": False}

    logger.info("Pipeline | %s | start", sample_id)

    qc_outdir = project_qc_dir(project_name) / sample_id
    qc_outdir.mkdir(parents=True, exist_ok=True)
    file_list = [r1] if r2 is None else [r1, r2]
    try:
        qc_out = run_qc_pipeline(sample_id, file_list, output_dir=str(qc_outdir))
        result["qc"] = qc_out
        logger.info("Pipeline | %s | QC OK", sample_id)
    except Exception as e:
        result["qc"] = {"error": str(e)}
        logger.error("Pipeline | %s | QC failed: %s", sample_id, e)

    qc_ok = result["qc"] and not result["qc"].get("error")

    if qc_ok and result["qc"].get("kraken2", {}).get("contamination"):
        _ng_pct = result["qc"]["kraken2"].get("ng_pct", 0.0)
        result["assembly"] = {
            "error": f"Excluded: only {_ng_pct:.1f}% of reads assigned to N. gonorrhoeae (< 80% threshold)",
        }
        result["excluded"] = True
        logger.warning("Pipeline | %s | excluded — ng_pct=%.1f%% < 80%%", sample_id, _ng_pct)
        return result

    r1_asm = Path(result["qc"]["trimmed_r1"]) if qc_ok and result["qc"].get("trimmed_r1") else r1
    _r2_trim = result["qc"].get("trimmed_r2") if qc_ok else None
    r2_asm = Path(_r2_trim) if _r2_trim and Path(_r2_trim).exists() else r2

    asm_outdir = project_assembly_dir(project_name) / sample_id
    try:
        asm = run_assembly_pipeline(
            sample_id=sample_id, r1=r1_asm, r2=r2_asm, outdir=asm_outdir,
            pre_trimmed=qc_ok,
        )
        result["assembly"] = {
            "contigs_path": str(asm.contigs) if asm.contigs else None,
            "logs":         asm.logs,
        }
        if asm.contigs:
            from backend.phylogeny.cgmlst import core_genome_gene_count
            try:
                core_genome_gene_count(str(asm.contigs))
            except Exception:
                logger.exception("Assembly | %s | core genome completeness check failed", sample_id)
            result["contigs_path"] = asm.contigs
            result["_contig_stats"] = _parse_contigs_stats(asm.contigs)
            logger.info("Pipeline | %s | assembly OK", sample_id)
        else:
            logger.error("Pipeline | %s | assembly produced no contigs", sample_id)
    except Exception as e:
        result["assembly"] = {"error": str(e)}
        logger.error("Pipeline | %s | assembly failed: %s", sample_id, e)

    if result["contigs_path"]:
        try:
            amr_out = run_amr_variant_calling(
                sample_id=sample_id,
                contigs=result["contigs_path"],
            )
            result["amr"] = amr_out
            logger.info("Pipeline | %s | AMR OK", sample_id)
        except Exception as e:
            result["amr"] = {"error": str(e)}
            logger.error("Pipeline | %s | AMR failed: %s", sample_id, e)

    return result


_FASTA_EXTS = {".fasta", ".fa", ".fna"}


def _process_single_fasta_sample(project_name: str, sample_id: str, fasta_path: Path) -> dict:
    assembled_path = fasta_path

    try:
        qc_out = run_qc_pipeline_fasta(
            sample_id, fasta_path, output_dir=str(project_qc_dir(project_name) / sample_id)
        )
        save_qc(project_name, sample_id, qc_out)
    except Exception as e:
        qc_out = {"error": str(e)}

    from backend.phylogeny.cgmlst import core_genome_gene_count
    try:
        core_genome_gene_count(str(fasta_path))
    except Exception:
        logger.exception("Assembly | %s | core genome completeness check failed", sample_id)
    stats = _parse_contigs_stats(fasta_path)
    stats_dict = stats.to_dict() if stats else {}
    save_assembly(project_name, sample_id, str(fasta_path), stats_dict)

    qc_ok = qc_out and not qc_out.get("error")
    if qc_ok and qc_out.get("kraken2", {}).get("contamination"):
        _ng_pct = qc_out["kraken2"].get("ng_pct", 0.0)
        logger.warning("Pipeline | %s | excluded (FASTA) — ng_pct=%.1f%% < 80%%", sample_id, _ng_pct)
        return {
            "qc":            qc_out,
            "assembly":      {"contigs_path": str(fasta_path), "logs": ""},
            "amr":           {"error": f"Excluded: only {_ng_pct:.1f}% of reads assigned to N. gonorrhoeae (< 80% threshold)"},
            "contigs_path":  assembled_path,
            "_contig_stats": stats,
            "excluded":      True,
        }

    try:
        amr_out = run_amr_variant_calling(sample_id=sample_id, contigs=fasta_path)
        save_amr(project_name, sample_id, amr_out)
    except Exception as e:
        amr_out = {"error": str(e)}

    return {
        "qc":              qc_out,
        "assembly":        {"contigs_path": str(fasta_path), "logs": ""},
        "amr":             amr_out,
        "contigs_path":    assembled_path,
        "_contig_stats":   stats,
        "excluded":        False,
    }


def run_full_pipeline(
    project_name: str,
    input_paths: list[Path],
    on_sample_done=None,
) -> dict:
    """
    Main entry point for the Full Pipeline: QC -> Assembly -> AMR ->
    Phylogeny for FASTQ samples. Pre-assembled FASTA inputs skip
    Assembly.

    Args:
        project_name: project to save results under.
        input_paths: uploaded FASTQ/FASTA files for the batch.
        on_sample_done: optional callback(sample_id, sample_result) fired
            as each sample finishes, for progress reporting.

    Returns:
        dict with one sub-dict per stage: qc, assembly, amr, phylogeny.
    """
    init_project(project_name)
    results: dict = {"qc": {}, "assembly": {}, "amr": {}, "phylogeny": {}}

    samples = _group_samples(input_paths)
    assembled_contigs: dict[str, Path] = {}

    fasta_samples = {
        sid: r for sid, r in samples.items()
        if r["R1"] is not None and Path(r["R1"]).suffix.lower() in _FASTA_EXTS
    }
    fastq_samples = {
        sid: r for sid, r in samples.items()
        if sid not in fasta_samples and r["R1"] is not None
    }

    if fasta_samples:
        _fasta_workers = max(1, min(len(fasta_samples), (os.cpu_count() or 2) // 2, 3))
        logger.info(f"FASTA pipeline executor: {_fasta_workers} workers")
        with ThreadPoolExecutor(max_workers=_fasta_workers) as fasta_executor:
            fasta_futures = {
                fasta_executor.submit(
                    _process_single_fasta_sample,
                    project_name,
                    sample_id,
                    reads["R1"],
                ): sample_id
                for sample_id, reads in fasta_samples.items()
            }
            for future in as_completed(fasta_futures):
                sample_id = fasta_futures[future]
                fasta_path = fasta_samples[sample_id]["R1"]
                try:
                    sr = future.result()
                except Exception as e:
                    logger.error(f"FASTA sample {sample_id} failed: {e}")
                    sr = {"qc": None, "assembly": None, "amr": {"error": str(e)}, "contigs_path": fasta_path, "excluded": False}
                results["amr"][sample_id] = sr["amr"]
                if sr.get("qc"):
                    results["qc"][sample_id] = sr["qc"]
                if sr.get("assembly"):
                    results["assembly"][sample_id] = sr["assembly"]
                if not sr.get("excluded"):
                    assembled_contigs[sample_id] = sr["contigs_path"]
                if on_sample_done is not None:
                    try:
                        on_sample_done(sample_id, sr)
                    except Exception:
                        pass

    try:
        _cg_limit    = int(open("/sys/fs/cgroup/memory.max").read().strip())
        _cg_current  = int(open("/sys/fs/cgroup/memory.current").read().strip())
        _headroom_gb = (_cg_limit - _cg_current) / 1e9
        _avail_gb    = _headroom_gb
        _ram_workers = max(1, int(_headroom_gb / 2.0))
        _ram_source  = f"cgroup headroom {_headroom_gb:.1f} GB / 2.0 GB per worker"
    except Exception:
        if _psutil is not None:
            _avail_gb    = _psutil.virtual_memory().available / 1e9
            _ram_workers = max(1, int(_avail_gb / 4.0))
            _ram_source  = f"psutil available {_avail_gb:.1f} GB / 4.0 GB per worker"
        else:
            _avail_gb    = 0.0
            _ram_workers = 1
            _ram_source  = "unknown (psutil unavailable)"
    _cpu_workers = max(1, (os.cpu_count() or 2) // 2)
    max_workers = min(len(fastq_samples), _ram_workers, _cpu_workers) if fastq_samples else 1
    logger.info("Pipeline executor: %d workers (%s)", max_workers, _ram_source)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _process_single_sample,
                project_name,
                sample_id,
                reads["R1"],
                reads["R2"],
            ): sample_id
            for sample_id, reads in fastq_samples.items()
        }
        for future in as_completed(futures):
            sample_id = futures[future]
            try:
                sr = future.result()
            except Exception as e:
                sr = {
                    "qc":          {"error": str(e)},
                    "assembly":    {"error": str(e)},
                    "amr":         {"error": str(e)},
                    "contigs_path": None,
                }

            if sr["qc"]:
                results["qc"][sample_id] = sr["qc"]
            if sr["assembly"]:
                results["assembly"][sample_id] = sr["assembly"]
            if sr["amr"]:
                results["amr"][sample_id] = sr["amr"]
            if sr["contigs_path"]:
                assembled_contigs[sample_id] = sr["contigs_path"]

            if sr["qc"] and not sr["qc"].get("error"):
                save_qc(project_name, sample_id, sr["qc"])
            if sr["contigs_path"]:
                _stats = sr.get("_contig_stats") or _parse_contigs_stats(sr["contigs_path"])
                save_assembly(project_name, sample_id, str(sr["contigs_path"]), _stats)
            if sr["amr"] and isinstance(sr["amr"], AMRResult):
                save_amr(project_name, sample_id, sr["amr"])

            if on_sample_done is not None:
                try:
                    on_sample_done(sample_id, sr)
                except Exception:
                    pass

    contig_list = list(assembled_contigs.values())
    if contig_list:
        try:
            phy = build_distance_tree([str(p) for p in contig_list])
            results["phylogeny"] = {
                "message":        phy["message"],
                "run_id":         phy["run_id"],
                "tree_path":      phy["tree_path"],
                "clusters":       phy.get("clusters", {}),
                "cluster_report": phy.get("cluster_report", {}),
                "matrix_path":    phy.get("matrix_path"),
                "upload_names":   phy.get("upload_names", []),
                "n_samples":      len(contig_list),
                "cgmlst_result":  phy.get("cgmlst_result"),
            }
            save_phylogeny(project_name, {
                "run_id":         phy["run_id"],
                "tree_path":      phy["tree_path"],
                "matrix_path":    phy.get("matrix_path"),
                "n_samples":      len(contig_list),
                "cluster_report": phy.get("cluster_report", {}),
                "cgmlst_result":  phy.get("cgmlst_result"),
            })
        except Exception as e:
            results["phylogeny"] = {"error": str(e), "tree_path": None}
    else:
        results["phylogeny"] = {"error": "No assemblies succeeded.", "tree_path": None}

    return results
