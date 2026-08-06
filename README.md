# NG-genomic-surveillance-platform
Integrated end-to-end pipeline for *N. gonorrhoeae* genomic surveillance, combining QC, assembly, AMR prediction and phylogenetic context with an interactive Streamlit interface.

## Overview

This project aims to provide a reproducible genomic surveillance platform for *Neisseria gonorrhoeae* that bridges raw sequencing data and public-health-ready outputs.

The platform is organised in independent modules:

1. Data cleaning & QC
2. De novo assembly & assembly QC
3. AMR determinant detection & genotype-to-phenotype mapping
4. Phylogenetic analysis & contextualisation
5. Interactive visualisation in Streamlit

## Pipeline modules


## Prepare environment (Conda)

conda env create -f env/environment.yml
conda activate ng_surv


## Folder structure


## Running the Streamlit app (local)

cd app/streamlit_app
streamlit run app.py

## Docker deployment (new machine)

```bash
cp docker-compose.dist.yml docker-compose.yml
cp .env.example .env
```

Edit `.env` to match the machine's available RAM. The worker memory limit controls how many samples are assembled in parallel — the pipeline reads the cgroup limit and divides by 2 GB per concurrent job:

| Machine RAM | `NG_WORKER_MEMORY` | Concurrent samples |
|-------------|--------------------|--------------------|
| 8 GB        | `4g`               | 1–2                |
| 16 GB       | `8g`               | 3–4                |
| 32 GB       | `16g`              | 7–8                |

```bash
docker compose up --build -d
```

The first build takes several minutes. Named volumes (`ng_projects`, `ng_results`, `ng_cgmlst_schema`, etc.) persist analysis results and downloaded schemas across rebuilds. To reset all data: `docker compose down -v`.
