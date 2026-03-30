# tb_mix -- Mixed infection detection for TB
# Converted from: containers/def/tb_mix.def
# Build context: project root (tb-lite/)

FROM python:3.12-slim

LABEL maintainer="TB-Lite pipeline"

ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    PATH=/opt:$PATH \
    PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive

COPY assets/scripts/tb_mix.py  /opt/tb_mix.py

RUN set -e \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
       build-essential procps coreutils gawk grep sed findutils \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir pandas pysam \
    && chmod +x /opt/tb_mix.py \
    && ln -s /opt/tb_mix.py /usr/local/bin/tb_mix.py
