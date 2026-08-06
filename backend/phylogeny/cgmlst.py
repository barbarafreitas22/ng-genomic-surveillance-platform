from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

try:
    import psutil as _psutil
except ImportError:
    _psutil = None

from ..paths import CGMLST_SCHEMA_DIR as SCHEMA_DIR

_PUBMLST_BASE  = "https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef"
_NG_CGMLST_V1  = 62

TRANSMISSION_THRESHOLD = 7
GENOGROUP_THRESHOLD    = 400

PROFILE_CACHE_DIR = SCHEMA_DIR / "profile_cache"


def is_chewbbaca_available() -> bool:
    return bool(shutil.which("chewBBACA.py") or shutil.which("chewbbaca"))


def is_schema_ready() -> bool:
    short = SCHEMA_DIR / "short"
    return short.is_dir() and any(short.glob("*.fasta"))


def _chewbbaca_bin() -> str:
    return shutil.which("chewBBACA.py") or shutil.which("chewbbaca") or "chewBBACA.py"


def _cpu_count() -> int:
    return max(os.cpu_count() or 1, 1)


def _available_memory_gb() -> float:
    try:
        _cg_limit   = int(open("/sys/fs/cgroup/memory.max").read().strip())
        _cg_current = int(open("/sys/fs/cgroup/memory.current").read().strip())
        return (_cg_limit - _cg_current) / 1e9
    except Exception:
        return (_psutil.virtual_memory().available / 1e9) if _psutil is not None else -1.0


def _chewbbaca_cpu_count(gb_per_cpu: float = 1.0) -> int:
    avail = _available_memory_gb()
    if avail < 0:
        return _cpu_count()
    return max(1, min(_cpu_count(), int(avail / gb_per_cpu)))


def _file_md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _pubmlst_get(path: str) -> dict:
    url = f"{_PUBMLST_BASE}/{path}"
    with urllib.request.urlopen(url, timeout=30) as r:
        return json.loads(r.read())


def _fetch_locus_fasta(args: tuple) -> tuple[str, str | None]:
    import time
    locus, dest = args
    fasta_path = Path(dest) / f"{locus}.fasta"
    if fasta_path.exists():
        return locus, None
    url = f"{_PUBMLST_BASE}/loci/{locus}/alleles/1"
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                d = json.loads(r.read())
            seq = d.get("sequence", "")
            if seq:
                fasta_path.write_text(f">{locus}_1\n{seq}\n")
                return locus, None
            return locus, "empty sequence"
        except urllib.error.HTTPError as exc:
            if exc.code == 429:
                time.sleep(2 ** attempt)
                continue
            return locus, str(exc)
        except Exception as exc:
            return locus, str(exc)
    return locus, "max retries (429)"


def download_schema() -> None:
    if not is_chewbbaca_available():
        raise RuntimeError(
            "chewBBACA not found in PATH. Install with:\n"
            "  conda install -c bioconda chewbbaca"
        )

    SCHEMA_DIR.mkdir(parents=True, exist_ok=True)
    loci_dir = SCHEMA_DIR / "loci_fasta"
    loci_dir.mkdir(exist_ok=True)

    data = _pubmlst_get(f"schemes/{_NG_CGMLST_V1}/loci")
    loci = [u.rstrip("/").split("/")[-1] for u in data.get("loci", [])]
    if not loci:
        raise RuntimeError("No loci returned from PubMLST for N. gonorrhoeae cgMLST v1.0.")

    args = [(loc, str(loci_dir)) for loc in loci]
    failed: list[str] = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futs = {pool.submit(_fetch_locus_fasta, a): a[0] for a in args}
        for fut in as_completed(futs):
            loc, err = fut.result()
            if err:
                failed.append(f"{loc}: {err}")

    present = len(list(loci_dir.glob("*.fasta")))
    if present < len(loci) * 0.95:
        raise RuntimeError(
            f"Only {present}/{len(loci)} loci downloaded successfully. "
            f"First failures: {failed[:5]}"
        )

    prep_out = SCHEMA_DIR / "prep_out"
    if prep_out.exists():
        shutil.rmtree(prep_out)

    proc = subprocess.run(
        [
            _chewbbaca_bin(), "PrepExternalSchema",
            "-g", str(loci_dir),
            "-o", str(prep_out),
            "--cpu", str(_chewbbaca_cpu_count()),
        ],
        capture_output=True, text=True, timeout=3600,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"chewBBACA PrepExternalSchema failed:\n{proc.stderr[-2000:]}"
        )

    if not (prep_out / "short").is_dir():
        raise RuntimeError(
            "PrepExternalSchema completed but short/ not found inside prep_out/."
        )

    for item in prep_out.iterdir():
        dst = SCHEMA_DIR / item.name
        if item.is_dir():
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(item, dst)
        else:
            shutil.copy2(item, dst)

    for hidden in prep_out.glob(".*"):
        shutil.copy2(hidden, SCHEMA_DIR / hidden.name)

    if not is_schema_ready():
        raise RuntimeError("Schema setup completed but not ready.")


def run_allele_calling(assembly_paths: list[str], output_dir: str) -> dict[str, dict[str, int]]:
    PROFILE_CACHE_DIR.mkdir(parents=True, exist_ok=True)

    cached_profiles: dict[str, dict[str, int]] = {}
    uncached_paths: list[str] = []
    path_to_md5: dict[str, str] = {}

    for p in assembly_paths:
        md5 = _file_md5(p)
        path_to_md5[p] = md5
        cache_file = PROFILE_CACHE_DIR / f"{md5}.json"
        if cache_file.exists():
            try:
                cached_profiles[Path(p).stem] = json.loads(cache_file.read_text())
                continue
            except Exception:
                pass
        uncached_paths.append(p)

    new_profiles: dict[str, dict[str, int]] = {}
    if uncached_paths:
        out = Path(output_dir) / "chewbbaca_out"
        out.mkdir(parents=True, exist_ok=True)

        list_file = out / "input_files.txt"
        list_file.write_text("\n".join(str(p) for p in uncached_paths) + "\n")

        cmd = [
            _chewbbaca_bin(), "AlleleCall",
            "-i", str(list_file),
            "-g", str(SCHEMA_DIR),
            "-o", str(out),
            "--cpu", str(_chewbbaca_cpu_count()),
            "--no-inferred",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if proc.returncode != 0:
            raise RuntimeError(f"chewBBACA AlleleCall failed:\n{proc.stderr[-2000:]}")

        results_tsv = next(out.rglob("results_alleles.tsv"), None)
        if results_tsv is None:
            raise RuntimeError(
                "Check the schema."
            )

        df = pd.read_csv(results_tsv, sep="\t", index_col=0)

        def _to_int(v):
            v = str(v).strip()
            if v.isdigit():
                return int(v)
            if v.startswith("INF-") and v[4:].isdigit():
                return 1_000_000 + int(v[4:])
            return 0

        for sample in df.index:
            sid = Path(str(sample)).stem
            profile = {locus: _to_int(df.loc[sample, locus]) for locus in df.columns}
            new_profiles[sid] = profile

            original_path = next((p for p in uncached_paths if Path(p).stem == sid), None)
            if original_path:
                md5 = path_to_md5.get(original_path)
                if md5:
                    try:
                        (PROFILE_CACHE_DIR / f"{md5}.json").write_text(json.dumps(profile))
                    except Exception:
                        pass

    return {**cached_profiles, **new_profiles}


def core_genome_completeness(assembly_path: str, output_dir: str) -> dict | None:
    if not (is_chewbbaca_available() and is_schema_ready()):
        return None
    try:
        profiles = run_allele_calling([assembly_path], output_dir)
    except Exception:
        return None
    profile = profiles.get(Path(assembly_path).stem)
    if not profile:
        return None
    total = len(profile)
    found = sum(1 for v in profile.values() if v != 0)
    return {
        "core_genes_found": found,
        "core_genes_total": total,
        "core_completeness_pct": round(found / total * 100, 1) if total else 0.0,
    }


def cached_core_genome_completeness(assembly_path: str) -> dict | None:
    cache_file = PROFILE_CACHE_DIR / f"{_file_md5(assembly_path)}.json"
    if not cache_file.exists():
        return None
    try:
        profile = json.loads(cache_file.read_text())
    except Exception:
        return None
    total = len(profile)
    found = sum(1 for v in profile.values() if v != 0)
    return {
        "core_genes_found": found,
        "core_genes_total": total,
        "core_completeness_pct": round(found / total * 100, 1) if total else 0.0,
    }


def _core_loci(profiles: dict, min_presence: float = 0.95) -> list[str]:
    if not profiles:
        return []
    all_loci = list(next(iter(profiles.values())).keys())
    if len(profiles) < 5:
        return all_loci
    n = len(profiles)
    return [
        locus for locus in all_loci
        if sum(1 for p in profiles.values() if p.get(locus, 0) != 0) / n >= min_presence
    ]


def compute_distance_matrix(
    profiles: dict[str, dict[str, int]],
    loci: list[str] | None = None,
) -> pd.DataFrame:
    samples = list(profiles.keys())
    if not samples:
        return pd.DataFrame()
    if loci is None:
        loci = list(profiles[samples[0]].keys())

    arr = np.array([[profiles[s].get(l, 0) for l in loci] for s in samples], dtype=np.int32)

    n = len(samples)
    mat = np.zeros((n, n), dtype=float)

    for i in range(n):
        va = arr[i]
        vb = arr[i + 1:]
        comparable = (va != 0) & (vb != 0)
        diffs = np.sum((va != vb) & comparable, axis=1).astype(float)
        mat[i, i + 1:] = diffs
        mat[i + 1:, i] = diffs

    return pd.DataFrame(mat, index=samples, columns=samples)


def cluster_at_threshold(dist_matrix: pd.DataFrame, threshold: int) -> dict[str, int]:
    samples = list(dist_matrix.index)
    if len(samples) == 0:
        return {}
    if len(samples) == 1:
        return {samples[0]: 1}

    condensed = squareform(dist_matrix.values)
    Z = linkage(condensed, method="single")
    labels = fcluster(Z, t=float(threshold), criterion="distance")
    return {s: int(lab) for s, lab in zip(samples, labels)}


def run_cgmlst_clustering(assembly_paths: list[str], output_dir: str) -> dict:
    profiles = run_allele_calling(assembly_paths, output_dir)
    if not profiles:
        raise RuntimeError("No allele profiles produced.")

    core = _core_loci(profiles)
    core_set = set(core)
    core_profiles = {
        sid: {l: v for l, v in p.items() if l in core_set}
        for sid, p in profiles.items()
    }

    dist = compute_distance_matrix(core_profiles, core)

    matrix_csv = str(Path(output_dir) / "cgmlst_distances.csv")
    dist.to_csv(matrix_csv)

    outbreak  = cluster_at_threshold(dist, TRANSMISSION_THRESHOLD)
    genogroup = cluster_at_threshold(dist, GENOGROUP_THRESHOLD)

    sample_stats: dict = {}
    for sid, alleles in core_profiles.items():
        assigned = sum(1 for v in alleles.values() if v != 0)
        pct      = round(assigned / len(core) * 100, 1) if core else 0.0
        sample_stats[sid] = {
            "loci_assigned":    assigned,
            "loci_total":       len(core),
            "pct_assigned":     pct,
            "outbreak_cluster": outbreak.get(sid),
            "genogroup":        genogroup.get(sid),
        }

    pairs: list[dict] = []
    samples = list(dist.index)
    for i, a in enumerate(samples):
        for b in samples[i + 1:]:
            d = int(dist.loc[a, b])
            pairs.append({
                "sample_a":       a,
                "sample_b":       b,
                "allele_diff":    d,
                "outbreak_link":  d <= TRANSMISSION_THRESHOLD,
                "same_genogroup": d <= GENOGROUP_THRESHOLD,
            })

    return {
        "core_loci_count":    len(core),
        "sample_stats":       sample_stats,
        "outbreak_clusters":  outbreak,
        "genogroup_clusters": genogroup,
        "pairwise":           pairs,
        "matrix_csv":         matrix_csv,
    }
