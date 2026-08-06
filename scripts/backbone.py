import os
import subprocess
import requests

OUTDIR    = "./data/phylogeny/genomes"
FASTQ_DIR = "./backend/pt_fastqc/downloads"
THREADS   = 4

os.makedirs(OUTDIR,    exist_ok=True)
os.makedirs(FASTQ_DIR, exist_ok=True)


# Samples and their resistance profiles 

SAMPLES = {
    "ERS4270824": "rCIP",
    "ERS4270828": "rCIP",
    "ERS4270893": "rAZM",
    "ERS4270984": "rAZM",
    "ERS4270985": "rCEF",
    "ERS4271005": "rPEN",
    "ERS4271062": "rPEN",
    "ERS4271063": "rTET",
    "ERR1560856": "rTET",
    "ERS4270815": "rSUL",
    "ERR1560803": "dPEN_dTET",
    "ERS4271240": "dTET",
    "ERS4271242": "dCEF",
    "ERS4271252": "dCEF",
    "ERS4271231": "FullySusceptible",
    "ERS4271237": "FullySusceptible",
    "ERS4271140": "rRI",
    "ERS4271133": "multiR",
    "ERS4271207": "multiR",
    "ERS4271239": "CEF_AZM",
    "ERS4270981": "dAZM",
    "ERS4271049": "dAZM_dCEF",
    "ERS4271178": "rSUL_dTET",
    "ERS4271289": "dTET",
}


def get_run_accession(accession):
    if accession.startswith("ERR"):
        return accession
    url = (
        f"https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={accession}&result=read_run"
        f"&fields=run_accession&format=tsv"
    )
    r = requests.get(url, timeout=30)
    lines = r.text.strip().splitlines()
    if len(lines) < 2:
        return None
    return lines[1].strip()


def run(cmd):
    print(f"  $ {' '.join(cmd)}")
    result = subprocess.run(cmd, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")



for accession, profile in SAMPLES.items():
    print(f"\n{'='*55}")
    print(f"  {accession}  ({profile})")
    print(f"{'='*55}")
    try:
        err = get_run_accession(accession)
    except Exception as e:
        print(f"  ERROR resolving accession: {e} — skipping")
        continue

    if not err:
        print(f"  ERROR: no run found for {accession} — skipping")
        continue

    print(f"  Run: {err}")

    sample_name = f"{err}_{profile}"
    assembly    = os.path.join(OUTDIR, f"{sample_name}.fna")

    if os.path.exists(assembly):
        print(f"  Already assembled — skipping")
        continue

    # Download fastq 
    fastq_out = os.path.join(FASTQ_DIR, err)
    os.makedirs(fastq_out, exist_ok=True)

    print(f"  Downloading {err}...")
    try:
        run(["fastq-dump", "--split-files", "--gzip", "--outdir", fastq_out, err])
    except RuntimeError as e:
        print(f"  ERROR downloading: {e} — skipping")
        continue

    r1 = os.path.join(fastq_out, f"{err}_1.fastq.gz")
    r2 = os.path.join(fastq_out, f"{err}_2.fastq.gz")

    if not os.path.exists(r1) or not os.path.exists(r2):
        print(f"  ERROR: fastq files not found — skipping")
        continue

    #  Assemble
    spades_out = os.path.join(FASTQ_DIR, f"{err}_assembly")
    print(f"  Assembling with SPAdes...")
    try:
        run([
            "spades.py",
            "--pe1-1", r1,
            "--pe1-2", r2,
            "-o", spades_out,
            "--threads", str(THREADS),
        ])
    except RuntimeError as e:
        print(f"  ERROR assembling: {e} — skipping")
        continue

    contigs = os.path.join(spades_out, "contigs.fasta")
    if not os.path.exists(contigs):
        print(f"  ERROR: contigs.fa not found — skipping")
        continue

    # Rename headers and save
    count = 0
    with open(contigs) as fin, open(assembly, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                count += 1
                fout.write(f">{sample_name}_contig{count}\n")
            else:
                fout.write(line)

    print(f"  Done: {assembly}  ({count} contigs)")

print(f"\n{'='*55}")
print(f"All done! Genomes saved to: {OUTDIR}")
print(f"{'='*55}\n")
