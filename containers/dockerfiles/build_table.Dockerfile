# build_table -- Python scripts for final table assembly
# Converted from: containers/def/build_table.def
# Build context: project root (tb-lite/)

FROM python:3.12-slim

LABEL maintainer="TB-Lite pipeline"

ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    PATH=/opt:$PATH \
    PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive

COPY assets/scripts/build_final_table.py  /opt/build_final_table.py
COPY assets/scripts/build_metrics_table.py  /opt/build_metrics_table.py
COPY assets/scripts/filter_tbmix.py  /opt/filter_tbmix.py
COPY assets/scripts/dr_parser.py  /opt/dr_parser.py
COPY assets/scripts/merge_ann_tables.py  /opt/merge_ann_tables.py

RUN set -e \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
       build-essential procps coreutils gawk grep sed findutils \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir pandas xlsxwriter openpyxl \
    && chmod +x /opt/build_final_table.py \
    && chmod +x /opt/build_metrics_table.py \
    && chmod +x /opt/filter_tbmix.py \
    && chmod +x /opt/dr_parser.py \
    && chmod +x /opt/merge_ann_tables.py \
    && ln -s /opt/build_final_table.py /usr/local/bin/build_final_table.py \
    && ln -s /opt/build_metrics_table.py /usr/local/bin/build_metrics_table.py \
    && ln -s /opt/filter_tbmix.py /usr/local/bin/filter_tbmix.py \
    && ln -s /opt/dr_parser.py /usr/local/bin/dr_parser.py \
    && ln -s /opt/merge_ann_tables.py /usr/local/bin/merge_ann_tables.py
