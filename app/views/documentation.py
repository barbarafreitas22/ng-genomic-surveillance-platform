from __future__ import annotations

import streamlit as st


def render() -> None:
    st.markdown("""
<div style="padding:0.6rem 0 1rem 0;">
  <div style="font-size:0.7rem;font-weight:700;text-transform:uppercase;letter-spacing:0.1em;color:#64748b;margin-bottom:0.25rem;">Documentation</div>
  <div style="font-size:1.6rem;font-weight:800;color:#0f172a;">Platform Technical Documentation</div>
</div>
""", unsafe_allow_html=True)
    st.markdown("""
    This page describes the design and implementation of each analysis module: tools, parameters,
    database sources, and construction procedures. It is intended as a reproducibility reference
    and methodological supplement.
    """)

    st.markdown("---")
    st.markdown("## Pipeline Architecture")
    st.markdown("""
    The platform implements a linear four-stage genomic analysis workflow:
    """)

    st.components.v1.html("""
    <div style="font-family:'Courier New',monospace; font-size:0.88rem;
                background:#f8fafc; color:#1e293b;
                border:1px solid #cbd5e1; border-radius:10px;
                padding:1.4rem 1.8rem; line-height:2.4;">
      <b>Raw FASTQ reads</b><br>
      &nbsp;&nbsp;&nbsp;&nbsp;↓&nbsp;&nbsp;fastp &nbsp;·&nbsp; Kraken2 &nbsp;·&nbsp; Coverage vs WHO&nbsp;F<br>
      <b>Module 1: QC &amp; Assembly</b><br>
      &nbsp;&nbsp;&nbsp;&nbsp;↓&nbsp;&nbsp;Trimmed reads<br>
      <b>Module 1 (cont.): De Novo Assembly</b>&nbsp;&nbsp;(SPAdes)<br>
      &nbsp;&nbsp;&nbsp;&nbsp;↓&nbsp;&nbsp;Assembled contigs (FASTA)<br>
      <b>Module 2: AMR Profiling</b>&nbsp;&nbsp;(Minimap2 · BCFtools · MLST · European 2020 rules)<br>
      &nbsp;&nbsp;&nbsp;&nbsp;↓&nbsp;&nbsp;Contigs<br>
      <b>Module 3: Phylogenetic Analysis</b>&nbsp;&nbsp;(SKA2 · RapidNJ)<br>
      &nbsp;&nbsp;&nbsp;&nbsp;↓<br>
      <b>Interactive phylogenetic tree + AMR report</b>
    </div>
    """, height=440)

    st.markdown("""
    Each module can also be run independently via the sidebar. All tools run inside a
    **Docker container** built from a single `Dockerfile` (micromamba), ensuring
    reproducibility across environments.

    | Component | Version / source |
    |---|---|
    | fastp | Bioconda |
    | SPAdes | Bioconda |
    | Minimap2 | Bioconda |
    | samtools / bcftools | Bioconda |
    | SKA2 (`ska2`) | Bioconda |
    | rapidNJ | Bioconda |
    | Kraken2 | Bioconda |
    | chewBBACA | Bioconda |
    | pyngoST | pip |
    | Python | 3.10 (conda-forge) |
    | Streamlit | pip |
    """)

    st.markdown("---")
    st.markdown("## Module 1: QC & Assembly")

    st.markdown("### fastp")
    st.markdown("""
    fastp performs adapter trimming and quality filtering on Illumina reads.
    Both single-end (SE) and paired-end (PE) modes are supported, detected
    automatically from the number of uploaded files.

    **Parameters:**

    | Flag | Value | Effect |
    |---|---|---|
    | `--thread` | 2 | CPU threads |
    | `--cut_right` | N/A | Sliding-window quality trimming from the 3′ end |
    | `--cut_window_size` | 4 | Window size for quality scan |
    | `--cut_mean_quality` | 20 | Clip window when mean Phred < 20 |
    | `--length_required` | 50 | Discard reads shorter than 50 bp after trimming |
    | `--detect_adapter_for_pe` | N/A | Auto-detect adapters for paired-end (PE only) |

    For paired-end runs, adapters are detected automatically without requiring an
    adapter file. For single-end runs, fastp uses its built-in adapter database.
    fastp also writes `fastp.json` and `fastp.html` quality reports alongside the
    trimmed reads. A standalone FastQC pass on the raw reads was dropped, its only
    two programmatically used fields, pre-trimming read count and GC%, are already
    present in fastp's own `before_filtering` summary, so running FastQC separately
    only duplicated numbers fastp already reports.
    """)

    st.markdown("### Species Confirmation with Kraken2")
    st.markdown("""
    Kraken2 is used to confirm that sequenced reads belong to *Neisseria gonorrhoeae*
    before assembly. A small custom database, built from *N. gonorrhoeae* strain diversity
    plus close *Neisseria* relatives, is used instead of the full RefSeq bacterial database
    to minimise memory footprint and classification time.
    """)

    with st.expander("Kraken2 database construction"):
        st.markdown("""
        **Build the database**
        ```bash
        python scripts/build_kraken2_db.py
        ```
        This downloads up to 100 complete *N. gonorrhoeae* genomes (plus the bundled WHO
        reference strains), up to 100 genomes each for five close relatives (*N. meningitidis*,
        *N. lactamica*, *N. cinerea*, *N. subflava*, *N. mucosa*), and up to 5 genomes each for
        six non-*Neisseria* outgroups common in genital/respiratory specimens (*Kingella kingae*,
        *Eikenella corrodens*, *Moraxella catarrhalis*, *Haemophilus influenzae*,
        *Staphylococcus aureus*, *Escherichia coli*), via the NCBI `datasets` CLI. `kraken2-build`
        then adds them to the library and builds the database, followed by `kraken2-build --clean`
        to drop the intermediate library FASTA and taxonomy download files. The close relatives
        and outgroups let off-target reads resolve to the correct genus/species instead of just
        "unclassified", and give contamination flags a concrete species label instead of none.

        The resulting database (~307 genomes, ~78 MB) consists of three binary files:
        - `hash.k2d`: k-mer hash table
        - `opts.k2d`: build options and parameters
        - `taxo.k2d`: taxonomy tree

        `data/kraken2_ng/` is small enough to bake directly into the Docker image
        (`COPY . .` in the `Dockerfile`, which runs this script during the build):
        no separate volume or per-machine download step is needed.

        **Classification command (per sample)**
        ```bash
        kraken2 --db /workspace/data/kraken2_ng \\
                --report kraken2_report.txt \\
                --output kraken2_output.txt \\
                --threads 4 \\
                input.fastq.gz
        ```
        Samples where *N. gonorrhoeae* does not appear as the dominant classification
        are flagged in the platform’s **Alerts** dashboard.
        """)

    st.markdown("---")
    st.markdown("## Module 1 (cont.): De Novo Assembly")

    st.markdown("### SPAdes")
    st.markdown("""
    SPAdes (St. Petersburg genome assembler) performs de novo assembly using a
    multi-k-mer de Bruijn graph approach. It is optimised for sequencing and 
    handles single-end and paired-end Illumina reads.

    **Parameters used:**

    | Parameter | Value | Notes |
    |---|---|---|
    | `--threads` | 4 | CPU threads |
    | `--memory` | dynamic | 35 % of available RAM, minimum 4 GB, maximum 32 GB |
    | Input mode | `--s1` (SE) or `--pe1-1` / `--pe1-2` (PE) | Detected automatically from uploaded file count |

    SPAdes automatically selects k-mer sizes based on read length. The assembler
    outputs `contigs.fasta` in the `spades_output/` subdirectory within the sample
    results folder.
    """)

    st.markdown("### Assembly Quality: Contiguity Metrics")
    st.markdown("""
    Contiguity statistics are computed directly from the assembly FASTA with Biopython:
    no external reference genome is used, since comparing a real clinical isolate against
    any single fixed reference strain would conflate genuine inter-strain divergence with
    assembly errors. The platform reports the following metrics directly in the UI:

    | Metric | Description |
    |---|---|
    | **Contigs** | Total number of assembled contigs |
    | **Total length** | Sum of all contig lengths (Mb) |
    | **N50 / N90** | Contig length at which 50% / 90% of the assembly is covered (kb) |
    | **L50 / L90** | Number of contigs needed to reach 50% / 90% of the assembly |
    | **auN** | Length-weighted mean contig length (Σ length² ÷ total length): a single-number contiguity summary less sensitive to assembly boundary effects than N50 alone |
    | **GC content** | GC percentage of the assembled genome |

    **Pass/caution/fail thresholds** follow the CDC AR Lab Network external quality
    assessment criteria for *N. gonorrhoeae* WGS [15]: at least 2 of the 3 metrics below
    must meet a tier's bar for the assembly to be scored at that tier.

    | Tier | Contigs | N50 | Total length |
    |---|---|---|---|
    | Pass | ≤150 | >30 kb | 2.0–2.2 Mb |
    | Caution | ≤180 | >20 kb | 1.8–2.2 Mb |

    Core genome completeness (BLASTN, below) and GC% (50–56%) are independent hard-fail
    checks specific to this platform, not part of the CDC criteria. The CDC source also
    defines a species-purity gate (≥90% pass / ≥85% caution reads identified as
    *N. gonorrhoeae*) and a coverage gate (≥40× pre-submission, ≥10× post-submission);
    this platform applies the ≥40× coverage bar during read QC (Module 1) and evaluates
    species purity separately via Kraken2, rather than merging them into this assembly-stage
    score. The CDC criteria also include an "MLST matches reference" gate, which applies
    only to CDC's own reference-isolate proficiency panel and has no equivalent for
    arbitrary clinical samples, so it is not implemented here.
    """)

    st.markdown("### Core Genome Completeness: BLASTN")
    st.markdown("""
    In addition to the contiguity metrics above, assembly completeness is assessed against a
    curated *N. gonorrhoeae* core gene set (1,713 genes), using BLASTN presence/absence
    screening (≥90% identity, ≥80% query coverage). The reported **Core genes %** is the
    percentage of these genes found present and intact in the assembly.

    This replaced an earlier genome-size-ratio metric (assembly length ÷ reference genome
    length), which could look deceptively high for an assembly that is the right *size* but
    missing or fragmenting real genes, and a later chewBBACA/cgMLST-based version, which reused
    the Module 3 clustering scheme (scheme 62, 1,638 loci) but added roughly two minutes of
    allele calling per genome to Module 1 QC. Samples below 90% core genome completeness are
    flagged in the assembly QC status, this threshold is checked upstream of AMR profiling,
    since a genome missing core loci may also be missing the resistance-panel genes needed for
    a reliable AMR call.
    """)

    st.markdown("---")
    st.markdown("## Module 2: AMR Profiling")

    st.markdown("### Resistance Gene Database")
    st.markdown("""
    Two sets of reference sequences were curated from NCBI RefSeq and the literature:
    """)

    st.markdown("**Chromosomal resistance genes**: variant calling (SNP detection)")
    st.markdown("""
    | Gene | Product | Resistance phenotype | Key mutations |
    |---|---|---|---|
    | `penA` | Penicillin-binding protein 2 (PBP2) | Ceftriaxone reduced susceptibility | A501V, A501P, A501T |
    | `gyrA` | DNA gyrase subunit A | Ciprofloxacin resistance | S91F, D95G, D95N, D95A |
    | `parC` | Topoisomerase IV subunit C | Ciprofloxacin resistance (with gyrA) | D86N, S87N, S87R |
    | `mtrR` | MtrCDE efflux pump repressor | Efflux pump overexpression | A39T |
    | `porB` | Outer membrane porin B | Reduced beta-lactam permeability | G120K, A121D, G120D |
    | `23SrRNA` | 23S ribosomal RNA | Azithromycin resistance | A2045G (high-level), C2597T (moderate) |
    | `blaTEM-1` | TEM-1 beta-lactamase | Penicillin resistance (chromosomal context) | N/A |
    """)

    st.markdown("**Plasmid resistance genes**. Presence/absence detection")
    st.markdown("""
    | Gene | Product | Resistance phenotype |
    |---|---|---|
    | `blaTEM-1` | TEM-1 beta-lactamase | Plasmid-mediated penicillin resistance (PPNG) |
    | `tet-M` | Tetracycline resistance protein M | Plasmid-mediated tetracycline resistance (TRNG) |

    Reference sequences were retrieved from NCBI GenBank using the `N. gonorrhoeae` strain FA1090
    (RefSeq AE004969) as the primary source for chromosomal genes, and from characterised
    resistance plasmids (pFAJ7 for *blaTEM-1*, pAV1 for *tet-M*).
    """)

    with st.expander("Gene database construction procedure"):
        st.markdown("""
        Each reference gene was extracted from the relevant NCBI entry and saved as a
        single-sequence FASTA file. FASTA index files (`.fai`) were generated with `samtools faidx`:

        ```bash
        # Example for penA
        samtools faidx data/genes/chromosomal/penA.fasta

        # Repeat for all gene sequences
        for f in data/genes/chromosomal/*.fasta data/genes/plasmid/*.fasta; do
            samtools faidx "$f"
        done
        ```

        The `.fai` index is required by BCFtools for variant calling against individual gene references.
        """)

    st.markdown("### Alignment: Minimap2")
    st.markdown("""
    Each assembled contig file is aligned against every reference gene sequence independently
    using Minimap2 with the `asm20` preset, which tolerates up to ~20% sequence divergence:
    needed for genes with distant or mosaic alleles (e.g. mosaic *penA*, mosaic *mtrR*/*mtrD*
    promoter blocks), where `asm5` (<5% divergence) would fail to align:

    ```bash
    minimap2 -a -x asm20 <ref_gene.fasta> <contigs.fasta>
    ```

    The SAM output is converted to sorted, indexed BAM with samtools:

    ```bash
    samtools view -bS alignment.sam | samtools sort -o sorted.bam
    samtools index sorted.bam
    ```
    """)

    st.markdown("### Variant Calling: BCFtools")
    st.markdown("""
    SNPs are called from the BAM alignment against each gene reference using a two-step
    BCFtools pipeline:

    ```bash
    bcftools mpileup -f <ref_gene.fasta> sorted.bam | bcftools call -mv -o variants.vcf
    ```

    - `mpileup`: piles up reads at each position to generate genotype likelihoods
    - `call -mv`: calls multiallelic SNPs and indels; `-v` outputs only variant sites

    Both SNPs and in-frame indels are parsed from the VCF. SNPs that fall in the same codon
    are combined before translation, so compound multi-nucleotide changes (e.g. three adjacent
    SNPs forming G120K) are reported as the correct single amino-acid substitution rather than
    three separate ones. In-frame insertions/deletions (length a multiple of 3) are reported as
    `ins<codon>` / `del<AA><codon>`; frameshift indels are discarded, as no clean codon-level
    notation applies.
    """)

    st.markdown("### Sequence Typing: MLST and NG-STAR")
    st.markdown("""
    Two independent PubMLST typing schemes are run against every assembly, using a shared
    BLAST-based allele-calling engine: each locus is BLASTed against its local allele-sequence
    set, and the best hit is classified as **exact** (≥99.9% identity, ≥95% coverage), **closest
    known** (≥95% identity, ≥80% coverage, reported with a `~` prefix), or **new**. A sequence
    type (ST) is only assigned when *all* loci in the scheme resolve to an exact match against
    the cached PubMLST profile table: a single inexact or novel locus is enough to leave the
    ST undetermined (`?`) rather than guess the closest one.

    - **MLST** (7 loci: `abcZ`, `adk`, `aroE`, `fumC`, `gdh`, `pdhC`, `pgm`): the standard
      *Neisseria* multilocus sequence typing scheme, used for general strain identification
      and population structure.
    - **NG-STAR** [14] (7 loci: `penA`, `mtrR`, `porB`, `ponA`, `gyrA`, `parC`, `23SrRNA`): an
      AMR-focused typing scheme covering the same genes tracked by the resistance-mutation
      panel above, enabling comparison against globally reported NG-STAR clonal complexes.

    Both schemes cache their PubMLST profile tables locally (30-day TTL) to avoid querying
    the PubMLST REST API on every run.
    """)

    st.markdown("### European 2020 (IUSTI) Interpretation Rules")
    st.markdown("""
    Called variants are translated from nucleotide position to amino acid notation
    (e.g. position 271 C→T in *gyrA* → **S91F**) using the standard genetic code.
    For RNA genes (`23SrRNA`), nucleotide notation is used directly (e.g. **A2045G**).

    Each detected mutation is matched against a curated rule table aligned with the
    **European (IUSTI/EuroSTI) 2020 Gonorrhoea Management Guidelines**:
    """)

    st.markdown("""
    | Mutation key | Phenotype | Therapy implication |
    |---|---|---|
    | penA A501V/P/T | Ceftriaxone reduced susceptibility | Avoid ceftriaxone/cefixime → gentamicin 240mg + azithromycin 2g |
    | 23SrRNA A2045G | High-level azithromycin resistance | Avoid azithromycin → ceftriaxone 1g IM monotherapy |
    | 23SrRNA C2597T | Moderate azithromycin resistance | Ceftriaxone 1g IM + azithromycin 2g (monitor MIC) |
    | gyrA S91F/D95G/N/A | Ciprofloxacin resistance | Avoid ciprofloxacin → ceftriaxone 1g IM + azithromycin 2g |
    | parC D86N/S87N/R | Ciprofloxacin resistance (additive) | Avoid ciprofloxacin |
    | mtrR A39T | Efflux pump overexpression | Monitor MIC |
    | porB G120K/A121D/G120D | Reduced beta-lactam permeability | Monitor MIC |
    | blaTEM-1 (plasmid) | PPNG: penicillin resistance | Avoid penicillin/ampicillin |
    | tet-M (plasmid) | TRNG: tetracycline resistance | Avoid tetracycline/doxycycline |
    """)

    st.markdown("### Recommended Treatment Regimens")
    st.markdown("""
    Treatment recommendations follow the **European (IUSTI/EuroSTI) 2020 Gonorrhoea Management
    Guidelines** (Unemo et al., 2020). Regimens are assigned per detected resistance phenotype
    and displayed in the clinical interpretation panel.

    **First-line: uncomplicated urogenital, anorectal infection (susceptibility unknown)** [1C]
    > Ceftriaxone **1g IM** (single dose) + azithromycin **2g orally** (single dose)

    | Clinical situation | Recommended regimen |
    |---|---|
    | Standard (susceptibility unknown) | Ceftriaxone 1g IM + azithromycin 2g orally |
    | Beta-lactam allergy | Spectinomycin 2g IM + azithromycin 2g orally [1B/1C] |
    | Beta-lactam allergy (alternative) | Gentamicin 240mg IM + azithromycin 2g orally [1B] |
    | IM injection contraindicated | Cefixime 400mg + azithromycin 2g orally [1B] |
    | Fluoroquinolone susceptibility confirmed | Ciprofloxacin 500mg orally (single dose) [1B] |
    | Ceftriaxone resistance | Gentamicin 240mg IM + azithromycin 2g orally [1B] |
    | Ceftriaxone resistance (alternative) | Spectinomycin 2g IM + azithromycin 2g orally [1B/1C] |

    **Notes:**
    - Azithromycin 2g should not be taken on an empty stomach. If GI side effects are anticipated,
      the dose may be split: azithromycin 1g + azithromycin 1g orally 6–12 h later.
    - Ciprofloxacin is contraindicated in pregnancy and should be used with caution in patients
      aged >60 years, with renal disease, or on corticosteroids.
    - Ceftriaxone 1g IM monotherapy is only an option where local ceftriaxone susceptibility
      data are comprehensive and recent, and TOC is mandatory.
    - Ceftriaxone resistance: three-site culture + AMR testing is required; notify public health
      authorities as mandated.
    - Co-infection with *C. trachomatis* is common in patients <30 years. If azithromycin is not
      used in the regimen, doxycycline 100mg bd for 7 days should be considered (not in pregnancy).
    """)

    st.markdown("""
    **Genomic Resistance Score**: a categorical score, not a sum of per-mutation weights.
    Detected mutations are first mapped to CDC phenotypes and, from these, to the affected
    antibiotic classes; the sample is then assigned to a single resistance category, and the
    category alone determines the score:

    | Category | Criterion | Score |
    |---|---|---|
    | Susceptible | No resistance phenotypes detected | 0.0 |
    | Low resistance | Only minor mechanisms (efflux upregulation, reduced susceptibility) | 0.2 |
    | Moderate resistance | ≥1 moderate-impact phenotype (e.g. ciprofloxacin resistance) | 0.4 |
    | High resistance | ≥1 high-impact phenotype (ceftriaxone or high-level azithromycin) | 0.6 |
    | MDR | ≥2 distinct antibiotic classes affected | 0.8 |
    | XDR | Ceftriaxone **and** azithromycin affected, plus ≥2 further classes | 1.0 |

    **This is an ordinal severity index, not a calibrated probability.** The category
    boundaries reflect real resistance biology (see below), but the even 0.2 spacing between
    them is an internal convention for ranking and sorting samples: it is not derived from any
    clinical outcome study, and the score does not mean "an X% chance of treatment failure."
    Only the ordering (susceptible < low < moderate < high < MDR < XDR) is meaningful; the
    numeric gaps between categories are not. Accordingly, the underlying function is named
    `resistance_severity_score()`, not a "probability" function.

    The MDR/XDR thresholds follow the MDR-GC / XDR-GC definitions of Unemo & Shafer (2014) [11]:
    XDR-GC requires resistance to both currently recommended therapies (ceftriaxone and
    azithromycin) plus ≥2 additional antibiotic classes: the combined ceftriaxone/high-level-
    azithromycin pattern reported in recent XDR case reports (Austria, France, Canada, Germany)
    [12] as the emerging threat of untreatable gonorrhoea. When counting affected classes,
    low-impact phenotypes (e.g. efflux-pump overexpression) are excluded, so that a single minor
    mechanism (which can trigger two related low-impact phenotypes at once) cannot by itself
    push a sample into the MDR category.

    This replaced an earlier per-mutation weighted-sum design: the individual numeric weights in
    that design had no direct literature source beyond the category a mutation already belonged
    to, so the score was simplified to use only the categorical structure, removing the
    uncited intermediate numbers. The score is used for risk stratification only; it does not
    replace clinical susceptibility testing.
    """)

    st.markdown("### Non-AMR Markers: Essential Gene Variants")
    st.markdown("""
    Shown as a separate, clearly-labelled informational section: these variants have **no
    effect** on the resistance category, phenotype calls, or the Genomic Resistance Score. They
    are detected with the same alignment and variant-calling engine as the resistance panel
    above (Minimap2 `asm20` + BCFtools), applied to nine genes with no known role in
    antimicrobial resistance, but with core functions in the bacterium:

    | Gene | Process | Product |
    |---|---|---|
    | `comA` | Natural transformation (competence) | DNA uptake across the outer membrane |
    | `ftsZ` | Cell division | Septum formation: tubulin homologue |
    | `recA` | DNA repair / recombination | Homologous recombination, SOS response |
    | `pilT`, `pilT2` | Motility (type IV pilus) | Twitching motility: retraction ATPase |
    | `tonB` | Iron acquisition | Energy transducer for TonB-dependent transporters |
    | `tbpB` | Iron acquisition | Transferrin-binding protein |
    | `fur` | Iron acquisition / stress response | Master regulator of iron uptake and oxidative stress |
    | `rpoH` | Stress response | Sigma-32 factor: activates heat-shock gene expression |

    Reference sequences were extracted directly from the *N. gonorrhoeae* FA1090 genome
    (GenBank AE004969.1), using GenBank annotation coordinates. This screen exists to
    characterise strain-level genomic diversity beyond the resistance phenotype: for example,
    whether variation in these genes correlates with the presence of resistance mutations,
    which would suggest broader selective pressure rather than neutral background variation.
    These are not clinically validated markers.
    """)

    st.markdown("---")
    st.markdown("## Module 3: Phylogenetic Analysis")

    st.markdown("### Reference Backbone")
    st.markdown("""
    To place uploaded genomes in a meaningful epidemiological context, a reference
    backbone of **33 genomes** was assembled from three sources:
    """)

    tab_who, tab_ena, tab_out = st.tabs(["WHO Reference Strains", "Clinical Isolates (ENA)", "Outgroup Genomes"])

    with tab_who:
        st.markdown("""
        **14 WHO reference strains** standardised panels for gonorrhoea antimicrobial
        susceptibility surveillance, used internationally for external quality assessment
        (WHO, 2016; Unemo et al., 2016):

        | Strain | Key resistance profile |
        |---|---|
        | WHO F | Chromosomally mediated resistance (CMRNG): penicillin |
        | WHO G | Fully susceptible |
        | WHO K | Penicillinase-producing (PPNG) |
        | WHO L | Fully susceptible |
        | WHO M | Penicillinase-producing (PPNG) |
        | WHO N | PPNG + tetracycline-resistant (TRNG) |
        | WHO O | Fully susceptible |
        | WHO P | Penicillinase-producing (PPNG) |
        | WHO U | PPNG + TRNG |
        | WHO V | Azithromycin reduced susceptibility |
        | WHO W | Ceftriaxone reduced susceptibility |
        | WHO X | Ceftriaxone + azithromycin reduced susceptibility |
        | WHO Y | Ciprofloxacin + azithromycin resistant |
        | WHO Z | Multi-drug resistant: ceftriaxone + azithromycin |

        Genome sequences retrieved from NCBI GenBank (BioProject PRJEB14414 and PRJNA342535).
        """)

    with tab_ena:
        st.markdown("""
        **21 clinical isolates** from the European Nucleotide Archive (ENA), selected to
        provide phenotypic representativeness across the major resistance profiles observed
        in European surveillance data. Isolates were assembled locally and annotated for
        AMR phenotype using the platform’s own AMR module:

        | ENA accession | AMR profile |
        |---|---|
        | ERR1560803 | Reduced penicillin + reduced tetracycline |
        | ERR1560856 | Tetracycline-resistant |
        | ERR3851425 | Reduced susceptibility (intermediate) |
        | ERR3851452 | Multi-drug resistant |
        | ERR3851478 | Sulfonamide-resistant |
        | ERR3851482 | Reduced azithromycin + ceftriaxone |
        | ERR3851551 | Fully susceptible |
        | ERR3851566 | Azithromycin-resistant |
        | ERR3851598 | Reduced tetracycline |
        | ERR3851600 | Ciprofloxacin-resistant |
        | ERR3851656 | Reduced ceftriaxone |
        | ERR3851665 | Ciprofloxacin-resistant |
        | ERR3851676 | Tetracycline-resistant |
        | ERR3851679 | Ceftriaxone-resistant |
        | ERR3851695 | Reduced azithromycin |
        | ERR3851796 | Penicillin-resistant |
        | ERR3851823 | Ceftriaxone + azithromycin resistant |
        | ERR3851871 | Reduced tetracycline |
        | ERR3851875 | Azithromycin-resistant |
        | ERR3851892 | Sulfonamide-resistant + reduced tetracycline |
        | ERR3851902 | Penicillin-resistant |

        Raw FASTQ files were downloaded from ENA using the `ena-fast-download` utility
        and assembled with SPAdes using the same parameters as Module 1.
        """)

    with tab_out:
        st.markdown("""
        **2 outgroup genomes**: phylogenetically close but distinct *Neisseria* species
        used to root the tree and provide evolutionary context:

        | Genome | Species | Source |
        |---|---|---|
        | `lactamica_genomic.fna` | *Neisseria lactamica* (CP031253) | NCBI RefSeq |
        | `mendinitis_genomic.fna` | *Neisseria meningitidis* (RM8376) | NCBI RefSeq |

        *N. lactamica* is a commensal Neisseria species with ~85% average nucleotide identity
        to *N. gonorrhoeae*, making it a suitable outgroup for intraspecific phylogenetic analyses.
        *N. meningitidis* serves as a more distant reference point.
        """)

    st.markdown("### SKA2: Split K-mer Analysis")
    st.markdown("""
    SKA2 (Split K-mer Analysis version 2) computes pairwise genomic distances between
    assemblies using split k-mers pairs of 15 bp flanking sequences with a central variable
    base. This approach is faster and more sensitive than alignment-based SNP methods for
    closely related isolates.

    **Workflow:**

    1. All genomes are concatenated into single-sequence FASTA files (multi-contig assemblies
       are linearised):
       ```python
       # All contigs joined into one sequence per sample
       combined_sequence = "".join(all_contig_sequences)
       ```

    2. A joint SKA2 sketch is built from all genomes in a single pass:
       ```bash
       ska build -o sketches -f file_list.tsv
       # file_list.tsv: tab-separated <sample_name> <fasta_path>
       ```

    3. All-vs-all pairwise distances are computed:
       ```bash
       ska distance sketches.skf > distances.tsv
       ```

    The distance matrix is converted to PHYLIP format for input to rapidNJ.

    The **pre-built backbone sketch** (`backbone.skf`) containing all reference and
    clinical isolate genomes is stored in `data/phylogeny/` and merged automatically
    with uploaded sample sketches on every run.
    """)

    st.markdown("### Tree Inference Method")
    st.markdown("""
    **Fast NJ: RapidNJ**
    Constructs a Neighbour-Joining tree directly from the SKA2 pairwise distance matrix.
    Suitable for local contextualisation within already-close groups and rapid outbreak screening.
    No bootstrap support values.
    ```bash
    rapidnj distances.phylip -i pd -o t > tree_nj.nwk
    ```
    """)

    st.markdown("### Genogroup Clustering")
    st.markdown("""
    **Genogroup (cgMLST, ≤400 allele differences).** Single-linkage clustering on the
    chewBBACA [13] cgMLST allele-difference matrix at a fixed threshold of **≤400 AD**,
    the same criterion used by Pathogenwatch for *N. gonorrhoeae* genomic clusters
    (PubMLST scheme 62, 1,638 loci). It is the only source of the tree's cluster coloring.

    **Methodological note.** Genogroup assignment requires chewBBACA, a locally available
    cgMLST schema, and an assembled genome. Module 3 only accepts assembled genomes
    (FASTA): raw FASTQ reads are uploaded in Module 1 (Quality Control & Assembly) or via
    the Full Pipeline, both of which assemble the reads before the genome reaches this
    clustering step, so every sample handed to Module 3 already qualifies for genogroup
    assignment. If chewBBACA or the schema is unavailable, the Genogroup field is left
    blank for that run. Because the threshold is absolute, cluster assignments are
    directly comparable across runs: adding or removing samples does not shift the
    cut-point.
    """)

    st.markdown("---")
    st.markdown("## Deployment with Docker")
    st.markdown("""
    The platform is orchestrated by `docker-compose.yml` and consists of two services:

    - **ng-platform**: Streamlit web interface, exposed on **port 8501** (6 GB memory
      limit). Serves the UI and coordinates job submission.
    - **ng-worker**: Background analysis worker that executes the bioinformatics pipeline.
      Uses the same Docker image as `ng-platform` but runs `python -m backend.worker`.
      Depends on `ng-platform` being healthy before starting.

    All bioinformatics tools are installed inside the Docker image via **micromamba**
    (from the `Dockerfile`). Reference databases in `data/` (gene FASTA files, Kraken2)
    are baked into the image at build time. `data/phylogeny/` is overridden at runtime
    with a bind mount, allowing updates to the backbone sketch, pyngoST, and cgMLST
    databases without rebuilding the image.

    **Volume mounts (ng-platform and ng-worker):**

    | Host / volume | Container path | Purpose |
    |---|---|---|
    | `./app` | `/workspace/app` | Application code (live-reload) |
    | `./backend` | `/workspace/backend` | Analysis modules (live-reload) |
    | `./data/phylogeny` | `/workspace/data/phylogeny` | Backbone sketch, pyngoST DB, cgMLST schema (all versioned in the repo) |
    | `ng_projects` (named) | `/workspace/app/projects` | Project data and sample uploads |
    | `ng_phylo_runs` (named) | `/workspace/backend/phylogeny/runs` | Phylogeny run outputs |
    | `ng_results` (named) | `/workspace/results` | Analysis result files |

    Named volumes (`ng_projects`, `ng_phylo_runs`, `ng_results`) are managed by Docker
    and persist across container recreations. The pyngoST database and cgMLST schema
    are not named volumes: they are ordinary files under `data/phylogeny/`, versioned
    in the repository, so a fresh `docker compose up --build` on a new machine already
    has them, no download step needed.
    """)

    with st.expander("Build and run commands"):
        st.markdown("""
        ```bash
        # Build image (first time or after Dockerfile changes)
        docker compose build

        # Start all services
        docker compose up -d

        # Rebuild without cache
        docker compose build --no-cache

        Access the platform
        # Platform:      http://localhost:8501
        ```
        """)

    st.markdown("---")
    st.markdown("## References")
    st.markdown("""
    [1] Bankevich A. et al. SPAdes: A New Genome Assembly Algorithm and Its Applications to
    Single-Cell Sequencing. *Journal of Computational Biology*, 19(5):455–477, 2012.
    doi: 10.1089/cmb.2012.0021

    [2] Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor.
    *Bioinformatics*, 34(17):i884–i890, 2018. doi: 10.1093/bioinformatics/bty560

    [3] Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*,
    34(18):3094–3100, 2018. doi: 10.1093/bioinformatics/bty191

    [4] Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2.
    *Genome Biology*, 20:257, 2019. doi: 10.1186/s13059-019-1891-0

    [5] Harris SR et al. SKA: Split k-mer analysis toolkit for bacterial genomic epidemiology.
    *Microbial Genomics*, 9(3), 2023. doi: 10.1099/mgen.0.000007

    [6] Simonsen M et al. Rapid neighbor-joining. *WABI 2008*, LNBI 5251:113–122, 2011.

    [7] Pinto M et al. Phylogenomics of *Neisseria gonorrhoeae* and the emergence and spread of
    antimicrobial resistance. *Microbial Genomics*, 6(7), 2020. doi: 10.1099/mgen.0.000481

    [8] Unemo M et al. The new WHO *Neisseria gonorrhoeae* reference strains for global quality
    assurance of laboratory investigations: phenotypic, genetic and reference genome
    characterization. *Journal of Antimicrobial Chemotherapy*, 71(11):3096–3108, 2016.
    doi: 10.1093/jac/dkw288

    [9] Unemo M et al. European guideline for the diagnosis and treatment of gonorrhoea in adults.
    *International Journal of STD & AIDS*, 31(12):1124–1134, 2020.
    doi: 10.1177/0956462420948791

    [10] Centers for Disease Control and Prevention. Sexually Transmitted Infections Treatment
    Guidelines, 2021: Updated 2024. *MMWR Recomm Rep*, 70(4):1–187, 2021.

    [11] Unemo M, Shafer WM. Antimicrobial resistance in *Neisseria gonorrhoeae* in the 21st
    century: past, evolution, and future. *Clinical Microbiology Reviews*, 27(3):587–613, 2014.
    doi: 10.1128/CMR.00010-14

    [12] Maubaret C, Caméléna F, Mrimèche M, et al. Two cases of extensively drug-resistant (XDR)
    *Neisseria gonorrhoeae* infection combining ceftriaxone-resistance and high-level azithromycin
    resistance, France, November 2022 and May 2023. *Eurosurveillance*, 28(37):2300456, 2023.
    doi: 10.2807/1560-7917.ES.2023.28.37.2300456

    [13] Silva M, Machado MP, Silva DN, et al. chewBBACA: A complete suite for gene-by-gene
    schema creation and strain identification. *Microbial Genomics*, 4(3), 2018.
    doi: 10.1099/mgen.0.000166

    [14] Demczuk W, Sidhu S, Unemo M, et al. *Neisseria gonorrhoeae* Sequence Typing for
    Antimicrobial Resistance (NG-STAR): a novel antimicrobial resistance multilocus typing
    scheme for tracking global dissemination of *N. gonorrhoeae* strains. *Journal of Clinical
    Microbiology*, 55(5), 2017. doi: 10.1128/JCM.00100-17: *volume/page details as recalled;
    verify against the published version before citing in the thesis.*

    [15] Reimche JL, Smith AC, Pham CD, Schmerer MW, Cartee JC, Bolden CB, Kersh EN, Gernert
    KM, et al. Establishing a *Neisseria gonorrhoeae* whole-genome sequencing external quality
    assessment for the Antimicrobial Resistance Laboratory Network (AR Lab Network) regional
    laboratories. *Microbiology Spectrum*, 14(3), 2026. doi: 10.1128/spectrum.02256-25
    """)


st.markdown("""

The system was implemented in the context of a Master’s thesis, with the aim of enabling rapid, automated, and standardized analysis of *N. gonorrhoeae* genomic data.
Source code and documentation are available on GitHub:
**https://github.com/barbarafreitas22/ng-genomic-surveillance-platform**
""")
