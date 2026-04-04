# rd -- runtime dependencies for region depth analysis
# The rd.py entrypoint is provided by Nextflow from the pipeline bin/ directory.

FROM python:3.11-slim

LABEL maintainer="TB-Lite pipeline"

ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive

RUN set -e \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
       build-essential procps \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir numpy pandas
