# NG-genomic-surveillance-platform

Integrated end-to-end pipeline for N. gonorrhoeae genomic surveillance, combining Quality Control, Assembly, AMR detection and phylogenetic context with an interactive Streamlit interface.

## Overview

This project aims to provide a reproducible genomic surveillance platform for *Neisseria gonorrhoeae* that bridges raw sequencing data and simplified outputs. The platform is organised into independent modules, each accessible from the Streamlit sidebar:

1. **Quality Control & Assembly**: read QC, *de novo* assembly, and assembly QC
2. **AMR Detection and Sequence Typing**: resistance determinant detection, genotype to phenotype mapping, and MLST / NG-STAR / NG-MAST typing
3. **Phylogenetic Analysis**: SNP and cgMLST based phylogenetic context
4. **Clinical Relevance**: *N. gonorrhoeae*-specific clinical background
5. **Technical Documentation**: methodology and tools behind each module

Per-sample epidemiological metadata (collection date, anatomical site, country, locality, health unit) can be added inline or imported from a CSV/TSV/Excel table (English or Portuguese headers) wherever samples are uploaded or analysed.

A **Full Pipeline Results** view runs all stages for a batch of samples and produces a shareable, self-contained HTML report plus a direct link to the project's results.

## Pipeline modules

| Stage | Tools |
|---|---|
| Read QC | fastp (trimming and read metrics), Kraken2 (species confirmation) |
| *De novo* assembly | SPAdes |
| Assembly QC | Biopython (contiguity metrics: N50, N90, L50, L90, auN, GC%), BLASTN (core genome completeness against a local panel of 1,713 *N. gonorrhoeae* core genes) |
| AMR profiling | Minimap2 + BCFtools variant calling against a curated resistance-gene panel; interpreted against European 2020 (IUSTI) treatment guidelines |
| Sequence typing | MLST, NG-STAR and NG-MAST: local BLAST against PubMLST allele sets, shared by the AMR and Phylogenetic modules |
| Phylogenetics | SKA2 (split k-mer, SNP distances), RapidNJ (tree construction), chewBBACA cgMLST (allele-based clustering) |

## Folder structure

```
.
├── app/                  Streamlit frontend
│   ├── app.py             entry point
│   └── views/              one module per page (qc_assembly, amr, phylogenetic, ...)
├── backend/               pipeline logic, decoupled from the UI
│   ├── qc.py, assembly.py, amr.py, mlst.py, ngstar.py, ngmast.py, ...
│   ├── worker.py           background job queue consumer
│   └── db.py                SQLite persistence per project
├── data/                  reference genome, gene panels, cgMLST/Kraken2 databases
├── config.ini             external tool binary paths (fastp, blastn, minimap2, ...)
├── scripts/                one-off setup scripts and admin utilities (schema downloads, job timing report)
├── tests/                  test fixtures
├── environment.yml        Conda environment definition
├── Dockerfile
└── docker-compose.yml      builds the image and runs the platform and the background worker
```

## Prepare environment (Conda)

```bash
conda env create -f environment.yml
conda activate ng
```

## Running the Streamlit app (local)

```bash
cd app
streamlit run app.py
```

## Docker deployment

Requirements: Git, [Git LFS](https://git-lfs.com) and Docker Desktop with at least 7 GB of memory
(Settings → Resources → Memory Limit).

```bash
git lfs install
git clone https://github.com/barbarafreitas22/ng-genomic-surveillance-platform.git
cd ng-genomic-surveillance-platform
git lfs pull
docker compose build --progress=plain
docker compose up -d
```

Then open http://localhost:8501.

- `git lfs install` must run before cloning: the Kraken2 database (`data/kraken2_ng/*.k2d`) is stored with Git LFS. Without it, these files are only small text pointers, the image builds without a working database, and species confirmation fails at runtime. `git lfs pull` fixes an existing clone; rebuild afterwards.
- The first build takes 15–30 minutes and needs internet access (conda packages, samtools source).
- If the build stops with `cannot allocate memory`, raise Docker Desktop's memory limit and run `docker compose build` again.

Useful commands:

```bash
docker compose ps          # container status
docker compose restart     # reload after code changes in app/ or backend/
docker compose down        # stop
docker compose down -v     # stop and delete all analysis results (named volumes)
```

`app/` and `backend/` are bind-mounted, so code changes only need `docker compose restart`; changes to the `Dockerfile` or `data/genes/` need a rebuild. Named volumes (`ng_projects`, `ng_results`, ...) keep analysis results across rebuilds.

## License

MIT: see [LICENSE](LICENSE).
