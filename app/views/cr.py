from __future__ import annotations

import streamlit as st


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
public health [1]. According to the World Health Organization (WHO), only in 2021 a total of 46,728 confirmed
cases of gonorrhea were reported in 27 countries of the European Union, reaffirming gonorrhea as the second
most prevalent sexually transmitted infection in Europe [2].

Beyond its clinical relevance, *N. gonorrhoeae* is characterized by a compact and highly adaptable genome,
consisting of a singular circular chromosome of approximately 2.2 megabase pairs (Mbp) encoding around
2,000 coding sequences [3].

---

## Microbiology

**Key biological features:**
- Gram negative diplococcus with characteristic kidney-shaped pairs
- Lacks a peptidoglycan layer thick enough for Gram-positive staining
- Highly adapted to evade host immunity through antigenic variation of surface structures
- Natural competence for DNA transformation facilitates rapid acquisition of resistance determinants

---

## Antimicrobial Resistance

The critical need for genomic surveillance is especially evident when examining the historical trajectory
of antimicrobial resistance (AMR) in *N. gonorrhoeae*. It has demonstrated a remarkable ability to develop
resistance to various classes of antibiotics used for its treatment, a consequence of the combination of
the bacterium's genetic adaptability and the indiscriminate use of multiple antibiotics to treat gonorrhea
[4, 5].

The complexity of mutations in *N. gonorrhoeae* makes their detection a demanding analytical task, largely
because most genetic determinants of AMR are located on the chromosome, and only the *blaTEM* gene [6] and
the *tetM* gene [7], which confer high level resistance, are known to be carried by plasmids in gonococci [8].
Resistance has been documented for:

| Antibiotic class | Agent | Mechanism |
|---|---|---|
| Beta-lactams | Penicillin | Plasmid-mediated (blaTEM-1) or chromosomal (penA, porB, mtrR) |
| Tetracyclines | Tetracycline / Doxycycline | Plasmid-mediated (tet-M) |
| Fluoroquinolones | Ciprofloxacin | Chromosomal (gyrA, parC) |
| Macrolides | Azithromycin | Chromosomal (23S rRNA A2045G, C2597T) |
| Extended-spectrum cephalosporins | Ceftriaxone | Chromosomal (penA mosaic alleles), last line therapy |

Currently, global gonorrhoea treatment guidelines recommend dual antimicrobial therapy: ceftriaxone, single
intramuscular dose, plus azithromycin, single oral dose, to mitigate resistance [2]. The European 2020
(IUSTI) guidelines followed by this platform recommend ceftriaxone 1g IM combined with azithromycin 2g
orally as the first line regimen, but reduced susceptibility to ceftriaxone remains an emerging global
concern. Monotherapy is recommended only in settings with robust surveillance confirming local
susceptibility; notably, there has been a global shift toward monotherapy in some regions, and since 2019
the UK and Japan have recommended a high dose of ceftriaxone alone [9] as a standard first line treatment
[10].

Given the imminent risk of therapeutic strategies becoming obsolete, and the fear that gonorrhoea will
evolve into an untreatable "superbug", specifically due to the emergence of multidrug-resistant (MDR) and
extensively drug-resistant (XDR) strains showing high level resistance to nearly all antimicrobial classes
[8], the WHO formally recognized this global threat, classifying the bacterium as a "Priority 2" (high
priority) microorganism [2]. In response to this urgency, the organization issued and updated guidelines,
in 2012, to raise awareness among healthcare professionals about the importance of clearly defining and
monitoring treatment failures [3].

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

[2] World Health Organization. WHO Guidelines for the Treatment of *Neisseria gonorrhoeae*. Geneva:
World Health Organization, 2016. ISBN 978-92-4-154971-4.
URL https://www.ncbi.nlm.nih.gov/books/NBK379221/.

[3] André Miguel Cascais Ferreira Pinto. Genome-scale approaches to strengthen *Neisseria gonorrhoeae*
epidemiological and antimicrobial resistance surveillance. PhD thesis. Supervised by João Paulo dos
Santos Gomes, Maria José Gonçalves Gaspar Borrego, and Luís Jaime Mota.

[4] Magnus Unemo and William M. Shafer. Antibiotic resistance in *Neisseria gonorrhoeae*: origin,
evolution, and lessons learned for the future. *Annals of the New York Academy of Sciences*,
1230:E19–E28, 2011. doi: 10.1111/j.1749-6632.2011.06215.x.

[5] Sarah J. Quillin and H. Steven Seifert. *Neisseria gonorrhoeae* host adaptation and pathogenesis.
*Nature Reviews Microbiology*, 16:226–240, 2018. doi: 10.1038/nrmicro.2017.169.

[6] W. A. Ashford, R. G. Golash, and V. G. Hemming. Penicillinase-producing *Neisseria gonorrhoeae*.
*The Lancet*, 308(7987):657–658, September 1976.

[7] Sandra A. Morse, Sylvia R. Johnson, John W. Biddle, and Marilyn C. Roberts. High level tetracycline
resistance in *Neisseria gonorrhoeae* is result of acquisition of streptococcal tetM determinant.
*Antimicrobial Agents and Chemotherapy*, 30(5):664–670, 1986.

[8] Magnus Unemo and William M. Shafer. Antimicrobial resistance in *Neisseria gonorrhoeae* in the 21st
century: past, evolution, and future. *Clinical Microbiology Reviews*, 27(3):587–613, 2014.
doi: 10.1128/CMR.00010-14.

[9] Magnus Unemo, Catriona S. Bradshaw, Jane S. Hocking, Henry J. C. de Vries, Suzanna C. Francis,
David Mabey, et al. Sexually transmitted infections: challenges ahead. *The Lancet Infectious Diseases*,
17(8):e235–e279, 2017. doi: 10.1016/S1473-3099(17)30310-9.

[10] Sancta St. Cyr, Lindley Barbee, Kimberly A. Workowski, et al. Update to CDC's Treatment Guidelines
for Gonococcal Infection, 2020. *MMWR Morbidity and Mortality Weekly Report*, 69(50):1911–1916, 2020.
doi: 10.15585/mmwr.mm6950a6.
    """)

