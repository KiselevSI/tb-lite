# ann_table -- Python scripts for annotation table processing
# Converted from: containers/def/ann_table.def
# Build context: project root (tb-lite/)

FROM python:3.12-slim

LABEL maintainer="TB-Lite pipeline"

ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    PATH=/opt:$PATH \
    PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive

# Copy scripts before installing deps (scripts change more often than deps)
COPY assets/scripts/add_name_strand.py  /opt/add_name_strand.py
COPY assets/scripts/dedup_ann_columns.py  /opt/dedup_ann_columns.py

RUN set -e \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
       build-essential procps coreutils gawk grep sed findutils \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir pandas xlsxwriter openpyxl \
    && chmod +x /opt/add_name_strand.py \
    && chmod +x /opt/dedup_ann_columns.py \
    && ln -s /opt/add_name_strand.py /usr/local/bin/add_name_strand.py \
    && ln -s /opt/dedup_ann_columns.py /usr/local/bin/dedup_ann_columns.py
