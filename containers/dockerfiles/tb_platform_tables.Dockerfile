# tb_platform_tables -- runtime dependencies for TB platform summary tables
# Table-generation scripts are provided by Nextflow from the pipeline bin/ directory.

FROM python:3.12-slim

LABEL maintainer="TB-Lite pipeline"

ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive

RUN set -e \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
       procps coreutils gawk grep sed findutils \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir pandas xlsxwriter openpyxl biopython
