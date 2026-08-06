from __future__ import annotations
import json
import os
import re
import tempfile
import time
import base64
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

import altair as alt
import pandas as pd
import streamlit as st

from views.shared import *


def render() -> None:
    st.markdown("""
<div style="padding:0.6rem 0 1rem 0;">
  <div style="font-size:0.7rem;font-weight:700;text-transform:uppercase;letter-spacing:0.1em;color:#64748b;margin-bottom:0.25rem;">Reference</div>
  <div style="font-size:1.6rem;font-weight:800;color:#0f172a;"><em>Neisseria gonorrhoeae</em> Clinical Relevance</div>
</div>
""", unsafe_allow_html=True)

    st.markdown("""
*Neisseria gonorrhoeae*, is a gram negative diplococcus, a strictly human pathogen and the causative agent
of gonorrhea, a sexually transmitted human infection that constitutes one of the greatest current threats to
public health [1].

---

## Microbiology

*N. gonorrhoeae* is a fastidious, obligate human pathogen that colonises mucosal surfaces of the
urogenital tract, rectum, pharynx, and conjunctiva. It is transmitted exclusively through direct
sexual contact and has no environmental reservoir.

**Key biological features:**
- Gram negative diplococcus with characteristic kidney-shaped pairs
- Lacks a peptidoglycan layer thick enough for Gram-positive staining
- Highly adapted to evade host immunity through antigenic variation of surface structures
- Natural competence for DNA transformation facilitates rapid acquisition of resistance determinants

---

## Antimicrobial Resistance

*N. gonorrhoeae* has demonstrated a remarkable capacity to acquire resistance to virtually every
antibiotic introduced for its treatment. Resistance has been documented for:

| Antibiotic class | Agent | Mechanism |
|---|---|---|
| Beta-lactams | Penicillin | Plasmid-mediated (blaTEM-1) or chromosomal (penA, porB, mtrR) |
| Tetracyclines | Tetracycline / Doxycycline | Plasmid-mediated (tet-M) |
| Fluoroquinolones | Ciprofloxacin | Chromosomal (gyrA, parC) |
| Macrolides | Azithromycin | Chromosomal (23S rRNA A2045G, C2597T) |
| Extended-spectrum cephalosporins | Ceftriaxone | Chromosomal (penA mosaic alleles) — last-line therapy |

The combination of the organism’s genetic adaptability and decades of sequential antibiotic exposure
has driven sequential resistance waves (Unemo & Shafer, 2011; Quillin & Seifert, 2018).
Ceftriaxone 1g IM combined with azithromycin 2g orally is the recommended first-line regimen (European 2020, IUSTI),
but reduced susceptibility to ceftriaxone is an emerging global concern.

---

## Epidemiology

Gonorrhoea is the second most commonly reported bacterial sexually transmitted infection globally.
In the EU/EEA, 106,331 cases were reported in 2024, a rate of 26.9 per 100,000 population and
an increase of over 300% since 2015 (ECDC Annual Epidemiological Report 2024).
The highest notification rates are recorded in Ireland, Malta, Iceland, Luxembourg, and Denmark.

---

## Role of Genomic Surveillance

Whole Genome sequencing (WGS) enables high resolution characterisation of *N. gonorrhoeae* beyond
what is achievable with conventional culture susceptibility testing. Genomic surveillance:

- Detects known resistance determinants and predicts phenotypic susceptibility
- Tracks transmission clusters and outbreak dynamics in near real time
- Identifies emerging resistant lineages before they become widespread
- Supports contact tracing and public health response at local and national levels
- Helps updating treatment guidelines locally and globally

This platform implements an automated WGS analysis workflow aligned with **European 2020 (IUSTI) guidelines**,
providing reproducible AMR profiling and phylogenetic analysis directly from raw sequencing data.

---

## References

[1] Benjamin I. Baarda, Ryszard A. Zielke, Alaina K. Holm, and Aleksandra E. Sikora. Comprehensive
Bioinformatic Assessments of the Variability of *Neisseria gonorrhoeae* Vaccine Candidates. *mSphere*,
6(1):e00977–20, February 2021. ISSN 2379-5042. doi: 10.1128/mSphere.00977-20.
URL https://journals.asm.org/doi/10.1128/mSphere.00977-20.
    """)

