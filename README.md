# NG-genomic-surveillance-platform

Integrated end-to-end pipeline for N. gonorrhoeae genomic surveillance, combining Quality Control, Assembly, AMR detection and phylogenetic context with an interactive Streamlit interface.

## Overview

This project aims to provide a reproducible genomic surveillance platform for *Neisseria gonorrhoeae* that bridges raw sequencing data and simplified outputs. The platform is organised into independent modules, each accessible from the Streamlit sidebar:

1. **Quality Control & Assembly**: read QC, *de novo* assembly, and assembly QC
2. **AMR Profiling**: resistance determinant detection and genotype to phenotype mapping
3. **Phylogenetic Analysis**: SNP and cgMLST based phylogenetic context
4. **Clinical Relevance**: *N. gonorrhoeae*-specific clinical background
5. **Technical Documentation**: methodology and tools behind each module

Per-sample epidemiological metadata (collection date, site, demographics) can be added inline wherever samples are uploaded or analysed.

A **Full Pipeline Results** view runs all stages for a batch of samples and produces a shareable, self-contained HTML report plus a direct link to the project's results.

## Pipeline modules

| Stage | Tools |
|---|---|
| Read QC | FastQC, fastp, Kraken2 (species confirmation) |
| *De novo* assembly | SPAdes |
| Assembly QC | Biopython (contiguity metrics: N50, N90, L50, L90, auN, GC%), BLASTN (core genome completeness against a local panel of 1,713 *N. gonorrhoeae* core genes) |
| AMR profiling | Minimap2 + BCFtools variant calling against a curated resistance-gene panel; interpreted against European 2020 (IUSTI) treatment guidelines |
| Sequence typing | MLST, NG-STAR and NG-MAST: local BLAST/minimap2 against PubMLST allele sets (AMR module); MLST and NG-STAR also available via pyngoST (Phylogenetic module) |
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
├── config.ini             external tool binary paths (fastqc, blastn, minimap2, ...)
├── scripts/                one-off setup scripts and admin utilities (schema downloads, job timing report)
├── tests/                  test fixtures
├── environment.yml        Conda environment definition
├── Dockerfile
├── docker-compose.yml      local dev compose file (builds the image)
└── docker-compose.dist.yml distributable compose file (currently builds locally; will switch to pulling a pre-built image once one is published)
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

## Docker deployment (new machine)

```bash
cp docker-compose.dist.yml docker-compose.yml
cp .env.example .env
```

Edit `.env` to match the machine's available RAM. `NG_WORKER_MEMORY` sets the worker container's hard memory limit; the worker reads this limit via cgroups and pauses claiming new samples whenever free headroom drops below 1.5 GB, so it self-throttles rather than crashing under memory pressure. As a practical starting point:


```bash
docker compose up --build -d
```

The first build takes several minutes. Named volumes (`ng_projects`, `ng_results`, etc.) persist analysis results across rebuilds. The cgMLST schema and pyngoST database are versioned under `data/phylogeny/` and bind-mounted, no download step needed on a fresh deploy. To reset non-versioned data: `docker compose down -v`.

## License

MIT: see [LICENSE](LICENSE).
