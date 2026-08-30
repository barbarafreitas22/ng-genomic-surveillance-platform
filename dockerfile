FROM mambaorg/micromamba:1.5.0

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

WORKDIR /workspace

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
        curl \
        default-jre-headless \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Split into multiple layers
RUN micromamba create -n ng -c conda-forge -c bioconda -c defaults \
        python=3.10 \
        pip \
    && micromamba clean --all --yes

RUN micromamba install -n ng -c conda-forge -c defaults \
        pandas numpy scipy scikit-learn biopython \
        streamlit altair matplotlib-base seaborn-base plotly \
        requests pyyaml python-dateutil rich psutil \
    && micromamba clean --all --yes

RUN micromamba install -n ng -c bioconda -c conda-forge \
        fastqc fastp \
        minimap2 samtools bcftools \
        kraken2 blast ncbi-datasets-cli \
        mlst mash rapidnj ska2=0.5.1 \
    && micromamba clean --all --yes

RUN micromamba install -n ng -c bioconda -c conda-forge --freeze-installed \
        spades \
    && micromamba clean --all --yes

# flash2 has no aarch64 conda package so it is built from source.
# samtools is built from source
# pins samtools to 0.1.x via perl-bio-samtools, preventing a conda-level upgrade.
# mlst uses BLAST for typing.
RUN apt-get update -qq && apt-get install -y --no-install-recommends \
        gcc make bzip2 zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev \
    && curl -fsSL https://github.com/dstreett/FLASH2/archive/refs/tags/2.2.00.tar.gz \
       | tar xz -C /tmp \
    && cd /tmp/FLASH2-2.2.00 && make && cp flash2 /opt/conda/envs/ng/bin/flash \
    && rm -rf /tmp/FLASH2-2.2.00 \
    && curl -fsSL https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 \
       | tar xj -C /tmp \
    && cd /tmp/samtools-1.21 \
    && ./configure --prefix=/opt/conda/envs/ng --without-curses --disable-plugins \
    && make -j4 && make install \
    && rm -rf /tmp/samtools-1.21 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN micromamba run -n ng pip install --no-cache-dir \
        vega-datasets newick pyfastx scikit-image pyngoST chewbbaca && \
    sed -i "s/from Bio.Align.Applications import MuscleCommandline/try:\n    from Bio.Align.Applications import MuscleCommandline\nexcept ImportError:\n    MuscleCommandline = None/" \
        /opt/conda/envs/ng/lib/python3.10/site-packages/pyngoST/pyngoST_utils.py

ENV PATH="/opt/conda/envs/ng/bin:${PATH}"
ENV PYTHONPATH="/workspace"

COPY . .

RUN true

RUN mkdir -p /workspace/.streamlit
COPY .streamlit/config.toml /workspace/.streamlit/config.toml

RUN mkdir -p \
    app/projects \
    app/uploads \
    app/results \
    backend/phylogeny/runs \
    data/phylogeny/pyngost_db \
    data/phylogeny/poppunk_db \
    data/phylogeny/cgmlst_schema \
    data/phylogeny/ska_cache/sample_sketch_cache \
    results

RUN micromamba run -n ng python scripts/build_kraken2_db.py \
    || echo "Kraken2 database build skipped — run scripts/build_kraken2_db.py manually after startup"

RUN chown -R mambauser:mambauser /workspace

USER mambauser

EXPOSE 8501

HEALTHCHECK --interval=30s --timeout=10s --retries=3 \
    CMD curl -f http://localhost:8501/_stcore/health || exit 1

ENTRYPOINT ["streamlit", "run", "app/app.py"]
