# ann_table -- runtime dependencies for SNP matrix table generation
# SnpEff database files are provided by Nextflow from params.snpeff_data_dir.

FROM mambaorg/micromamba:2.0.5

LABEL maintainer="TB-Lite pipeline"

USER root

ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/bin:$PATH

RUN set -e \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
       procps coreutils gawk grep sed findutils \
       default-jre-headless \
    && rm -rf /var/lib/apt/lists/*

RUN micromamba install -y -n base -c conda-forge -c bioconda \
       python=3.12 \
       pandas \
       bcftools \
       snpsift \
    && micromamba clean --all --yes

RUN printf '%s\n' \
       '#!/bin/sh' \
       'JAR=$(find /opt/conda -name SnpSift.jar -print -quit)' \
       '[ -n "$JAR" ] || { echo "SnpSift.jar not found under /opt/conda" >&2; exit 127; }' \
       'exec java -jar "$JAR" "$@"' \
    > /usr/local/bin/snpSift \
    && chmod +x /usr/local/bin/snpSift \
    && ln -s /usr/local/bin/snpSift /usr/local/bin/SnpSift

USER $MAMBA_USER
