# rd -- Region depth analysis script
# Converted from: containers/def/rd.def
# Build context: project root (tb-lite/)

FROM python:3.11-slim

LABEL maintainer="TB-Lite pipeline"

ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    PATH=/opt/rdscanner:$PATH \
    PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive

COPY assets/scripts/rd.py  /opt/rd.py

RUN set -e \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
       build-essential procps \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir numpy pandas \
    && chmod +x /opt/rd.py \
    && ln -s /opt/rd.py /usr/local/bin/rd.py
