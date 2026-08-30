import hashlib
import json
import logging
import os
import glob
import shutil
import subprocess
import tempfile
import threading
import uuid
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
from Bio import SeqIO

try:
    import psutil as _psutil
except ImportError:
    _psutil = None

from backend.phylogeny.run_pyngost import run_typing, is_pyngost_ready
from backend.qc import run_fastp as _run_fastp, FASTP_BIN as _FASTP_BIN

logger = logging.getLogger(__name__)

_SKETCH_LOCK = threading.Lock()


def _memory_aware_workers(max_by_count: int, gb_per_worker: float) -> int:
    """Cap a worker pool by both CPU count (via max_by_count) and available RAM."""
    try:
        _cg_limit   = int(open("/sys/fs/cgroup/memory.max").read().strip())
        _cg_current = int(open("/sys/fs/cgroup/memory.current").read().strip())
        _avail_gb   = (_cg_limit - _cg_current) / 1e9
    except Exception:
        _avail_gb = (_psutil.virtual_memory().available / 1e9) if _psutil is not None else -1.0
    if _avail_gb < 0:
        return 1
    return max(1, min(max_by_count, int(_avail_gb / gb_per_worker)))

from ..paths import (
    BASE_GENOMES_DIR, RUNS_DIR,
    BACKBONE_SKETCH, BACKBONE_HASH_FILE, BACKBONE_METADATA_PATH, BACKBONE_FILE_LIST,
    SAMPLE_SKETCH_CACHE_DIR, SAMPLE_SKETCH_MANIFEST,
)

_CPU_THREADS = str(min(os.cpu_count() or 4, 8))

BACKBONE_SKETCH        = str(BACKBONE_SKETCH)
BACKBONE_HASH_FILE     = str(BACKBONE_HASH_FILE)
BACKBONE_METADATA_PATH = str(BACKBONE_METADATA_PATH)
SAMPLE_SKETCH_CACHE_DIR = str(SAMPLE_SKETCH_CACHE_DIR)
SAMPLE_SKETCH_MANIFEST  = str(SAMPLE_SKETCH_MANIFEST)


def _load_backbone_metadata() -> dict:
    try:
        with open(BACKBONE_METADATA_PATH) as f:
            data = json.load(f)
        return {k: v for k, v in data.items() if not k.startswith("_")}
    except Exception:
        return {}


def ensure_dirs():
    os.makedirs(BASE_GENOMES_DIR, exist_ok=True)
    os.makedirs(RUNS_DIR, exist_ok=True)
    os.makedirs(SAMPLE_SKETCH_CACHE_DIR, exist_ok=True)


def run(cmd):
    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed:\n{' '.join(cmd)}\n\nSTDERR:\n{proc.stderr}"
        )
    return proc.stdout


def concatenate_fasta(input_path, output_path):
    sample_name = os.path.splitext(os.path.basename(input_path))[0]
    combined = "".join(str(record.seq) for record in SeqIO.parse(str(input_path), "fasta"))

    with open(output_path, "w") as out:
        out.write(f">{sample_name}\n{combined}\n")

    return sample_name


def get_base_genomes():
    exts = ["*.fa", "*.fasta", "*.fna", "*.fas"]
    files = []
    for ext in exts:
        files.extend(glob.glob(os.path.join(BASE_GENOMES_DIR, ext)))
    return sorted(os.path.abspath(f) for f in files)


def _backbone_fingerprint(genome_paths: list) -> str:
    h = hashlib.md5()
    for p in sorted(genome_paths):
        h.update(p.encode())
        try:
            h.update(str(os.path.getmtime(p)).encode())
        except OSError:
            pass
    return h.hexdigest()


def _get_cached_backbone_hash() -> str:
    try:
        with open(BACKBONE_HASH_FILE) as f:
            return f.read().strip()
    except FileNotFoundError:
        return ""


def _save_backbone_hash(fingerprint: str) -> None:
    with open(BACKBONE_HASH_FILE, "w") as f:
        f.write(fingerprint)


def _ensure_backbone_sketch(backbone_genomes: list, concat_dir: str) -> None:
    """Build backbone.skf only when the reference genomes have changed."""
    with _SKETCH_LOCK:
        fingerprint = _backbone_fingerprint(backbone_genomes)
        if (
            fingerprint == _get_cached_backbone_hash()
            and os.path.exists(BACKBONE_SKETCH)
        ):
            return

        backbone_list = str(BACKBONE_FILE_LIST)
        file_list_rows = []

        for genome in backbone_genomes:
            stem = os.path.splitext(os.path.basename(genome))[0]
            concat_path = os.path.join(concat_dir, f"{stem}.fna")
            concatenate_fasta(genome, concat_path)
            file_list_rows.append((stem, concat_path))

        with open(backbone_list, "w") as f:
            for name, path in file_list_rows:
                f.write(f"{name}\t{path}\n")

        prefix = os.path.splitext(BACKBONE_SKETCH)[0]
        run(["ska", "build", "-o", prefix, "-f", backbone_list, "--threads", _CPU_THREADS])
        _save_backbone_hash(fingerprint)


def _file_md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_sketch_manifest() -> dict:
    try:
        with open(SAMPLE_SKETCH_MANIFEST) as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def _save_sketch_manifest(manifest: dict) -> None:
    with open(SAMPLE_SKETCH_MANIFEST, "w") as f:
        json.dump(manifest, f, indent=2)


def _cached_sample_sketch(genome_path: str, node_name: str, concat_path: str, manifest: dict) -> str:
    """Return a cached per-sample .skf, rebuilding only when content or label changes."""
    abs_path = os.path.abspath(genome_path)
    current_md5 = _file_md5(genome_path)

    cached = manifest.get(abs_path, {})
    if (
        cached.get("md5") == current_md5
        and cached.get("node_name") == node_name
        and os.path.exists(cached.get("sketch_path", ""))
    ):
        return cached["sketch_path"]

    cache_key = hashlib.md5(f"{current_md5}:{node_name}".encode()).hexdigest()
    sketch_prefix = os.path.join(SAMPLE_SKETCH_CACHE_DIR, cache_key)
    sketch_path = sketch_prefix + ".skf"

    unique_prefix = f"{sketch_prefix}_{uuid.uuid4().hex}"
    tmp_list = f"{unique_prefix}_fl.tsv"
    with open(tmp_list, "w") as fh:
        fh.write(f"{node_name}\t{concat_path}\n")
    try:
        run(["ska", "build", "-o", unique_prefix, "-f", tmp_list, "--threads", "4"])
    finally:
        if os.path.exists(tmp_list):
            os.unlink(tmp_list)

    unique_sketch = unique_prefix + ".skf"
    if not os.path.exists(unique_sketch):
        raise RuntimeError(f"SKA build failed for {os.path.basename(genome_path)}")

    os.replace(unique_sketch, sketch_path)

    manifest[abs_path] = {"md5": current_md5, "node_name": node_name, "sketch_path": sketch_path}
    return sketch_path


def _cached_sample_sketch_fastq(r1: str, r2: str, node_name: str, manifest: dict, trim: bool) -> str:
    r1_abs, r2_abs = os.path.abspath(r1), os.path.abspath(r2)
    current_md5 = hashlib.md5((_file_md5(r1) + _file_md5(r2)).encode()).hexdigest()
    manifest_key = f"{r1_abs}|{r2_abs}"

    cached = manifest.get(manifest_key, {})
    if (
        cached.get("md5") == current_md5
        and cached.get("node_name") == node_name
        and cached.get("trim") == trim
        and os.path.exists(cached.get("sketch_path", ""))
    ):
        return cached["sketch_path"]

    cache_key = hashlib.md5(f"{current_md5}:{node_name}:{trim}".encode()).hexdigest()
    sketch_prefix = os.path.join(SAMPLE_SKETCH_CACHE_DIR, cache_key)
    sketch_path = sketch_prefix + ".skf"

    trim_dir = None
    if trim:
        trim_dir = tempfile.mkdtemp(dir=SAMPLE_SKETCH_CACHE_DIR, prefix=f"trim_{cache_key}_")
        r1_use, r2_use = _run_fastp(Path(r1), Path(r2), Path(trim_dir), _FASTP_BIN)
    else:
        r1_use, r2_use = r1, r2

    unique_prefix = f"{sketch_prefix}_{uuid.uuid4().hex}"
    tmp_list = f"{unique_prefix}_fl.tsv"
    with open(tmp_list, "w") as fh:
        fh.write(f"{node_name}\t{r1_use}\t{r2_use}\n")
    try:
        run(["ska", "build", "-o", unique_prefix, "-f", tmp_list, "--threads", "4"])
    finally:
        if os.path.exists(tmp_list):
            os.unlink(tmp_list)
        if trim_dir:
            shutil.rmtree(trim_dir, ignore_errors=True)

    unique_sketch = unique_prefix + ".skf"
    if not os.path.exists(unique_sketch):
        raise RuntimeError(f"SKA build failed for {node_name} (FASTQ)")

    os.replace(unique_sketch, sketch_path)

    manifest[manifest_key] = {
        "md5": current_md5, "node_name": node_name, "trim": trim, "sketch_path": sketch_path,
    }
    return sketch_path


def ska_pairwise_to_matrix(tsv_path):
    rows = []
    with open(tsv_path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            try:
                snps = float(parts[2])
                prop = float(parts[3]) if len(parts) > 3 else 0.0
            except ValueError:
                continue
            rows.append((parts[0], parts[1], snps, prop))

    samples = sorted(set([r[0] for r in rows] + [r[1] for r in rows]))
    snp_matrix  = pd.DataFrame(0.0, index=samples, columns=samples)
    prop_matrix = pd.DataFrame(0.0, index=samples, columns=samples)

    for a, b, snps, prop in rows:
        snp_matrix.loc[a, b]  = snps
        snp_matrix.loc[b, a]  = snps
        prop_matrix.loc[a, b] = prop
        prop_matrix.loc[b, a] = prop

    return snp_matrix, prop_matrix


def write_phylip(matrix, out_path):
    with open(out_path, "w") as f:
        f.write(f"{len(matrix)}\n")
        for s in matrix.index:
            row = " ".join(map(str, matrix.loc[s]))
            f.write(f"{s} {row}\n")


def build_distance_tree(uploaded_samples):
    """
    Main entry point for Module 3 and the Full Pipeline. Both callers only
    ever pass assembled genomes (FASTA): Module 3's own uploader is
    FASTA-only, and the Full Pipeline always assembles reads before
    handing genomes to this function. Sketches the genomes with SKA2,
    merges them with the pre-built backbone sketch, computes pairwise
    SNP distances, builds an NJ tree (RapidNJ). If chewBBACA/schema are
    available, layers on cgMLST genogroup clustering (see
    backend/phylogeny/cgmlst.py), the only source for `clusters`.

    Args:
        uploaded_samples: list of sample dicts or paths (assembled
            FASTA) to place on the tree.

    Returns:
        dict with tree_path, matrix_path, cluster_report (per-sample
        nearest-neighbour + typing + cgmlst stats), clusters (genogroup
        assignment per sample, empty for samples without one), and
        cgmlst_result (None if cgMLST wasn't available for this run).
    """
    ensure_dirs()

    run_id  = str(uuid.uuid4())
    run_dir = os.path.join(RUNS_DIR, run_id)
    os.makedirs(run_dir, exist_ok=True)

    uploads_dir = os.path.join(run_dir, "uploads")
    concat_dir  = os.path.join(run_dir, "concat")
    os.makedirs(uploads_dir, exist_ok=True)
    os.makedirs(concat_dir,  exist_ok=True)

    dist_tsv    = os.path.join(run_dir, "distances.tsv")
    dist_phylip = os.path.join(run_dir, "distances.phylip")
    treefile    = os.path.join(run_dir, "tree.nwk")

    backbone_genomes = get_base_genomes()
    _ensure_backbone_sketch(backbone_genomes, concat_dir)

    backbone_names = [
        os.path.splitext(os.path.basename(g))[0] for g in backbone_genomes
    ]

    sketch_manifest = _load_sketch_manifest()

    import re as _re
    entries: list = []
    seen_names: set = set()
    for idx, sample in enumerate(uploaded_samples):
        if isinstance(sample, (str, Path)):
            sample = {"kind": "fasta", "path": str(sample)}

        is_fasta = sample["kind"] == "fasta"
        if is_fasta:
            gpath = Path(sample["path"])
            base = sample.get("name") or gpath.stem
        else:
            base = sample["name"]

        node_name = _re.sub(r'[(),:;\s]+', '_', base).strip('_') or f"sample-{idx+1}"
        if node_name in seen_names:
            node_name = f"{node_name}_{idx + 1}"
        seen_names.add(node_name)

        if is_fasta:
            raw_dest = os.path.join(uploads_dir, node_name + gpath.suffix)
            shutil.copy2(sample["path"], raw_dest)
            concat_path = os.path.join(concat_dir, f"{node_name}.fna")
            concatenate_fasta(raw_dest, concat_path)
            entries.append({
                "kind": "fasta", "name": node_name,
                "concat_path": concat_path, "raw_path": raw_dest,
            })
        else:
            entries.append({
                "kind": "fastq", "name": node_name,
                "r1": sample["r1"], "r2": sample["r2"], "trim": sample.get("trim", True),
            })

    def _build_one(entry):
        if entry["kind"] == "fasta":
            return _cached_sample_sketch(entry["raw_path"], entry["name"], entry["concat_path"], sketch_manifest)
        return _cached_sample_sketch_fastq(entry["r1"], entry["r2"], entry["name"], sketch_manifest, entry["trim"])

    _cpu_cap = max(1, min(len(entries), (os.cpu_count() or 4) // 4))
    _sketch_workers = _memory_aware_workers(_cpu_cap, gb_per_worker=0.5)
    with ThreadPoolExecutor(max_workers=_sketch_workers) as _sketch_pool:
        upload_sketches = list(_sketch_pool.map(_build_one, entries))

    _save_sketch_manifest(sketch_manifest)

    upload_names = {e["name"] for e in entries}
    all_names = backbone_names + list(upload_names)
    if len(all_names) < 2:
        raise ValueError("At least 2 genomes are required.")

    if not entries:
        dist_output   = run(["ska", "distance", BACKBONE_SKETCH])
        merged_sketch = BACKBONE_SKETCH
    else:
        merged_prefix = os.path.join(run_dir, "merged")
        run(["ska", "merge", BACKBONE_SKETCH] + upload_sketches + ["-o", merged_prefix])
        merged_sketch = f"{merged_prefix}.skf"
        if not os.path.exists(merged_sketch):
            raise RuntimeError("SKA merge failed — merged.skf not created.")
        dist_output = run(["ska", "distance", merged_sketch])

    with open(dist_tsv, "w") as f:
        f.write(dist_output)

    matrix, prop_matrix = ska_pairwise_to_matrix(dist_tsv)
    write_phylip(matrix, dist_phylip)

    newick_out = run(["rapidnj", dist_phylip, "-i", "pd", "-o", "t"])
    if not newick_out.strip():
        raise RuntimeError("Empty tree returned by rapidNJ")
    with open(treefile, "w") as f:
        f.write(newick_out.strip() + "\n")

    upload_raw_paths = [e["raw_path"] for e in entries if e["kind"] == "fasta"]
    input_types = {
        e["name"]: ("assembly" if e["kind"] == "fasta"
                     else f"reads ({'trimmed' if e['trim'] else 'raw'})")
        for e in entries
    }

    clusters: dict = {}

    matrix_csv = os.path.join(run_dir, "distances.csv")
    matrix.to_csv(matrix_csv)

    backbone_meta = _load_backbone_metadata()
    cluster_report: dict = {}
    for sample in matrix.index:
        row    = matrix.loc[sample]
        others = row.drop(sample, errors="ignore")
        if others.empty:
            continue
        nn_name = others.idxmin()
        nn_dist = float(others.min())
        nn_pct  = round(float(prop_matrix.loc[sample, nn_name]) * 100, 4)
        entry = {
            "is_uploaded":       sample in upload_names,
            "input_type":        input_types.get(sample, "assembly"),
            "nearest_neighbour": nn_name,
            "nn_distance":       round(nn_dist, 1),
            "nn_distance_pct":   nn_pct,
            "nn_is_uploaded":    nn_name in upload_names,
        }
        if sample in backbone_meta:
            entry["backbone_info"] = backbone_meta[sample]
        if nn_name in backbone_meta:
            entry["nn_backbone_info"] = backbone_meta[nn_name]
        cluster_report[sample] = entry

    pyngost_ready = is_pyngost_ready()
    if upload_raw_paths and pyngost_ready:
        try:
            typing = run_typing(upload_raw_paths, run_dir)
            for sample, tdata in typing.items():
                if sample in cluster_report:
                    cluster_report[sample]["typing"] = tdata
        except Exception:
            logger.exception("MLST/NG-STAR typing failed for run_dir=%s", run_dir)

    cgmlst_result = None
    try:
        from backend.phylogeny.cgmlst import (
            is_chewbbaca_available, is_schema_ready, run_cgmlst_clustering,
        )
        if upload_raw_paths and is_chewbbaca_available() and is_schema_ready():
            cgmlst_result = run_cgmlst_clustering(upload_raw_paths, run_dir)
            cg_geo = cgmlst_result.get("genogroup_clusters", {})
            for sample, stats in cgmlst_result.get("sample_stats", {}).items():
                if sample in cluster_report:
                    cluster_report[sample]["cgmlst"] = stats
            clusters.update(cg_geo)
    except Exception:
        logger.exception("cgMLST clustering failed for run_dir=%s", run_dir)

    n_total = len(matrix)
    return {
        "message":        f"Tree built successfully with {n_total} genomes",
        "run_id":         run_id,
        "tree_path":      treefile,
        "clusters":       clusters,
        "cluster_report": cluster_report,
        "matrix_path":    matrix_csv,
        "upload_names":   list(upload_names),
        "pyngost_ready":  pyngost_ready,
        "cgmlst_result":  cgmlst_result,
    }
