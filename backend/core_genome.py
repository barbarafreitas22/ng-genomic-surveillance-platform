from __future__ import annotations

import configparser as _configparser
import hashlib
import json
import re
import subprocess

from .paths import CONFIG_PATH, CORE_GENOME_FASTA, RESULTS_DIR

_config = _configparser.ConfigParser()
_config.read(CONFIG_PATH)
BLASTN = _config["PATHS"].get("BLASTN", "blastn")

CORE_GENOME_CACHE_DIR = RESULTS_DIR / "core_genome_cache"

_CORE_MIN_PID = 90.0
_CORE_MIN_COV = 0.80


def _file_md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _core_genome_gene_total() -> int:
    if not CORE_GENOME_FASTA.exists():
        return 0
    with open(CORE_GENOME_FASTA) as f:
        return sum(1 for line in f if line.startswith(">"))


def _core_genome_genes() -> list[dict]:
    genes = []
    if not CORE_GENOME_FASTA.exists():
        return genes
    with open(CORE_GENOME_FASTA) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            tags = dict(re.findall(r"\[(\w+)=([^\]]*)\]", header))
            genes.append({
                "id":        header.split()[0].removeprefix("lcl|"),
                "gene":      tags.get("gene", ""),
                "locus_tag": tags.get("locus_tag", ""),
                "product":   tags.get("protein", ""),
            })
    return genes


def core_genome_gene_count(
    assembly_path: str, min_pid: float = _CORE_MIN_PID, min_cov: float = _CORE_MIN_COV,
    need_genes: bool = False,
) -> dict | None:
    """
    Alternative completeness metric against data/genes/core_genome.fasta Uses BLASTN
    presence/absence by identity+coverage instead of allele calling,
    since there's no multi-allele catalog to call against.

    Returns:
        dict with core_genes_found/total/core_completeness_pct, or None
        if the reference file is missing.
    """
    if not CORE_GENOME_FASTA.exists():
        return None

    total = _core_genome_gene_total()
    if not total:
        return None

    CORE_GENOME_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_file = CORE_GENOME_CACHE_DIR / f"{_file_md5(assembly_path)}.json"
    if cache_file.exists():
        try:
            cached = json.loads(cache_file.read_text())
            if not need_genes or "genes" in cached:
                return cached
        except Exception:
            pass

    proc = subprocess.run(
        [
            BLASTN,
            "-query",         str(CORE_GENOME_FASTA),
            "-subject",       str(assembly_path),
            "-outfmt",        "6 qseqid pident length qlen",
            "-perc_identity", str(min_pid),
            "-dust",          "no",
        ],
        capture_output=True, text=True, timeout=600,
    )
    if proc.returncode != 0:
        return None

    best_cov: dict[str, float] = {}
    best_pid: dict[str, float] = {}
    for line in proc.stdout.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        qid, pident, length, qlen = parts[0].removeprefix("lcl|"), parts[1], parts[2], parts[3]
        try:
            cov = int(length) / int(qlen)
        except (ValueError, ZeroDivisionError):
            continue
        if cov > best_cov.get(qid, 0.0):
            best_cov[qid] = cov
            best_pid[qid] = float(pident)

    found = sum(1 for c in best_cov.values() if c >= min_cov)
    result = {
        "core_genes_found": found,
        "core_genes_total": total,
        "core_completeness_pct": round(found / total * 100, 1) if total else 0.0,
        "genes": {
            qid: [round(best_pid[qid], 1), round(cov * 100, 1)]
            for qid, cov in best_cov.items()
        },
    }
    try:
        cache_file.write_text(json.dumps(result))
    except Exception:
        pass
    return result


def cached_core_genome_gene_count(assembly_path: str) -> dict | None:
    """Read-only lookup for core_genome_gene_count's cache, does not run BLASTN."""
    cache_file = CORE_GENOME_CACHE_DIR / f"{_file_md5(assembly_path)}.json"
    if not cache_file.exists():
        return None
    try:
        return json.loads(cache_file.read_text())
    except Exception:
        return None


def core_genome_gene_table(assembly_path: str, compute: bool = False) -> list[dict] | None:
    core = (
        core_genome_gene_count(assembly_path, need_genes=True) if compute
        else cached_core_genome_gene_count(assembly_path)
    )
    if not core or "genes" not in core:
        return None

    hits = core["genes"]
    rows = []
    for g in _core_genome_genes():
        pid, cov = hits.get(g["id"], (None, None))
        if cov is None:
            status = "Absent"
        elif cov >= _CORE_MIN_COV * 100:
            status = "Present"
        else:
            status = "Partial"
        rows.append({
            "Gene":         g["gene"] or "-",
            "Locus tag":    g["locus_tag"],
            "Product":      g["product"],
            "Status":       status,
            "Identity (%)": pid,
            "Coverage (%)": cov,
        })
    return rows
