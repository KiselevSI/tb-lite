# tblg -- TB lineage classification tool
# Converted from: containers/def/tblg.def
# Build context: project root (tb-lite/)

FROM python:3.11-slim

LABEL maintainer="TB-Lite pipeline"

ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    PATH=/opt/rdscanner:$PATH \
    PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive

RUN set -e \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
       build-essential procps \
    && rm -rf /var/lib/apt/lists/* \
    && pip install "numpy<2.0" "pandas<2.0" tblg
