# tb_mix -- runtime dependencies for mixed infection detection
# The tb_mix.py entrypoint is provided by Nextflow from the pipeline bin/ directory.

FROM python:3.12-slim

LABEL maintainer="TB-Lite pipeline"

ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive

RUN set -e \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
       build-essential procps coreutils gawk grep sed findutils \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir pandas pysam
