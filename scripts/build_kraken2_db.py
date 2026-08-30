import glob
import json
import logging
import os
import shutil
import subprocess
import sys
import zipfile

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s",
                     datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)

DB_DIR = "./data/kraken2_ng"
TMP_DIR = os.path.join(DB_DIR, "_build_tmp")
WHO_GENOMES_DIR = "./data/phylogeny/genomes"

TARGET_TAXON = "Neisseria gonorrhoeae"
TARGET_TAXID = "485"
TARGET_LIMIT = 100
RELATIVE_LIMIT = 100
RELATIVE_TAXA = [
    "Neisseria meningitidis",
    "Neisseria lactamica",
    "Neisseria cinerea",
    "Neisseria subflava",
    "Neisseria mucosa",
]

# Non-Neisseria outgroups, common co-isolates/contaminants in genital and
# oral/respiratory clinical specimens. Give the classifier something concrete
# to assign contamination reads to, instead of them falling to "unclassified".
OUTGROUP_LIMIT = 5
OUTGROUP_TAXA = [
    "Kingella kingae",
    "Eikenella corrodens",
    "Moraxella catarrhalis",
    "Haemophilus influenzae",
    "Staphylococcus aureus",
    "Escherichia coli",
]


def run(cmd):
    result = subprocess.run(cmd, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")


def download_limited_genomes(taxon, outdir, limit):
    """
    The `datasets` CLI no longer supports --limit (removed upstream), so the
    genome list is fetched first via `summary` and truncated client-side.
    """
    summary = subprocess.run(
        ["datasets", "summary", "genome", "taxon", taxon,
         "--assembly-level", "complete", "--as-json-lines"],
        capture_output=True, text=True,
    )
    if summary.returncode != 0:
        raise RuntimeError(f"datasets summary failed: {summary.stderr}")

    accessions = []
    for line in summary.stdout.strip().splitlines():
        try:
            accessions.append(json.loads(line)["accession"])
        except Exception:
            continue
        if len(accessions) >= limit:
            break

    if not accessions:
        raise RuntimeError(f"No complete genomes found for {taxon}")

    zip_path = os.path.join(TMP_DIR, f"{taxon.replace(' ', '_')}.zip")
    run([
        "datasets", "download", "genome", "accession", *accessions,
        "--include", "genome", "--filename", zip_path,
    ])
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(outdir)


def add_who_reference_strains(ng_dir):
    who_files = sorted(glob.glob(os.path.join(WHO_GENOMES_DIR, "WHO*_genomic.fna")))
    added = 0
    for src in who_files:
        label = os.path.basename(src).replace("_genomic.fna", "")
        dst = os.path.join(ng_dir, f"{label}.fna")
        with open(src) as fin, open(dst, "w") as fout:
            for line in fin:
                if line.startswith(">"):
                    fout.write(f">{label}|kraken:taxid|{TARGET_TAXID}\n")
                else:
                    fout.write(line)
        added += 1
    return added


def main():
    if os.path.exists(os.path.join(DB_DIR, "hash.k2d")):
        return

    os.makedirs(TMP_DIR, exist_ok=True)

    ng_dir = os.path.join(TMP_DIR, "ng")
    download_limited_genomes(TARGET_TAXON, ng_dir, TARGET_LIMIT)
    add_who_reference_strains(ng_dir)

    for taxon in RELATIVE_TAXA:
        rel_dir = os.path.join(TMP_DIR, taxon.replace(" ", "_"))
        download_limited_genomes(taxon, rel_dir, RELATIVE_LIMIT)

    for taxon in OUTGROUP_TAXA:
        out_dir = os.path.join(TMP_DIR, taxon.replace(" ", "_"))
        try:
            download_limited_genomes(taxon, out_dir, OUTGROUP_LIMIT)
        except RuntimeError as e:
            logger.warning("Skipping outgroup %s: %s", taxon, e)

    run(["kraken2-build", "--download-taxonomy", "--use-ftp", "--db", DB_DIR])

    added = 0
    for root, _dirs, files in os.walk(TMP_DIR):
        for f in files:
            if f.endswith(".fna"):
                run(["kraken2-build", "--add-to-library", os.path.join(root, f), "--db", DB_DIR])
                added += 1

    if added == 0:
        logger.error("No genomes were downloaded")
        sys.exit(1)

    run(["kraken2-build", "--build", "--db", DB_DIR, "--threads", str(os.cpu_count() or 4)])

    run(["kraken2-build", "--clean", "--db", DB_DIR])
    shutil.rmtree(TMP_DIR, ignore_errors=True)

    expected = ["hash.k2d", "opts.k2d", "taxo.k2d"]
    missing = [f for f in expected if not os.path.exists(os.path.join(DB_DIR, f))]
    if missing:
        logger.error("Build did not produce %s in %s", missing, DB_DIR)
        sys.exit(1)


if __name__ == "__main__":
    main()
