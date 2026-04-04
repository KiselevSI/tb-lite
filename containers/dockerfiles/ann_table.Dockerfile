# ann_table -- runtime dependencies for SNP matrix table generation
# Post-processing scripts are provided by Nextflow from the pipeline bin/ directory.

FROM python:3.12-slim

LABEL maintainer="TB-Lite pipeline"

ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive \
    SNPEFF_HOME=/opt/snpEff

COPY snpEff_latest_core/snpEff /opt/snpEff

RUN set -e \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
       default-jre-headless \
       bcftools \
       procps coreutils gawk grep sed findutils \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir pandas \
    && printf '#!/bin/sh\nexec java -jar /opt/snpEff/SnpSift.jar "$@"\n' > /usr/local/bin/snpSift \
    && chmod +x /usr/local/bin/snpSift
