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


TARGET_TAXON = "Neisseria gonorrhoeae"
TARGET_LIMIT = 15
RELATIVE_TAXA = [
    "Neisseria meningitidis",
    "Neisseria lactamica",
    "Neisseria cinerea",
    "Neisseria subflava",
]


def run(cmd):
    logger.info("$ %s", " ".join(cmd))
    result = subprocess.run(cmd, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")


def download_genomes(taxon, outdir, extra_args=None):
    zip_path = os.path.join(TMP_DIR, f"{taxon.replace(' ', '_')}.zip")
    cmd = ["datasets", "download", "genome", "taxon", taxon, "--include", "genome",
           "--filename", zip_path]
    cmd += extra_args or []
    run(cmd)
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(outdir)


def main():
    if os.path.exists(os.path.join(DB_DIR, "hash.k2d")):
        logger.info("Database already present at %s — skipping.", DB_DIR)
        return

    os.makedirs(TMP_DIR, exist_ok=True)

    ng_dir = os.path.join(TMP_DIR, "ng")
    download_genomes(TARGET_TAXON, ng_dir,
                      ["--assembly-level", "complete", "--limit", str(TARGET_LIMIT)])

    for taxon in RELATIVE_TAXA:
        logger.info("Downloading reference genome for %s...", taxon)
        rel_dir = os.path.join(TMP_DIR, taxon.replace(" ", "_"))
        download_genomes(taxon, rel_dir, ["--reference"])

    run(["kraken2-build", "--download-taxonomy", "--db", DB_DIR])

    added = 0
    for root, _dirs, files in os.walk(TMP_DIR):
        for f in files:
            if f.endswith(".fna"):
                run(["kraken2-build", "--add-to-library", os.path.join(root, f), "--db", DB_DIR])
                added += 1

    if added == 0:
        logger.error("No genomes were downloaded")
        sys.exit(1)

    logger.info("Building database from %d genomes...", added)
    run(["kraken2-build", "--build", "--db", DB_DIR, "--threads", str(os.cpu_count() or 4)])

    logger.info("Cleaning up intermediate library files...")
    run(["kraken2-build", "--clean", "--db", DB_DIR])
    shutil.rmtree(TMP_DIR, ignore_errors=True)

    expected = ["hash.k2d", "opts.k2d", "taxo.k2d"]
    missing = [f for f in expected if not os.path.exists(os.path.join(DB_DIR, f))]
    if missing:
        logger.error("Build did not produce %s in %s", missing, DB_DIR)
        sys.exit(1)

    logger.info("Done. Kraken2 database ready at %s", DB_DIR)


if __name__ == "__main__":
    main()
