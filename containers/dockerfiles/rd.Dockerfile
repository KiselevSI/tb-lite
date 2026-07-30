# rd -- runtime for RD (region of difference) scanning.
# The rd_scan.py entrypoint is provided by Nextflow from the pipeline bin/ directory
# and depends only on the Python standard library.

FROM python:3.11-slim

LABEL maintainer="TB-Lite pipeline"

ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive

RUN set -e \
    && apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
