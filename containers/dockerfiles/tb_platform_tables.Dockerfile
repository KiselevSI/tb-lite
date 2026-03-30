# tb_platform_tables -- Python scripts for TB platform summary tables
# Build context: project root (tb-lite/)

FROM python:3.12-slim

LABEL maintainer="TB-Lite pipeline"

ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    PATH=/opt:$PATH \
    PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive

COPY assets/scripts/filter_tbmix.py        /opt/filter_tbmix.py
COPY assets/scripts/build_general_table.py /opt/build_general_table.py
COPY assets/scripts/deletions_to_csv.py    /opt/deletions_to_csv.py
COPY assets/scripts/profiler_parser.py     /opt/profiler_parser.py

RUN set -e \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
       procps coreutils gawk grep sed findutils \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir pandas xlsxwriter openpyxl biopython \
    && chmod +x /opt/filter_tbmix.py \
    && chmod +x /opt/build_general_table.py \
    && chmod +x /opt/deletions_to_csv.py \
    && chmod +x /opt/profiler_parser.py \
    && ln -s /opt/filter_tbmix.py /usr/local/bin/filter_tbmix.py \
    && ln -s /opt/build_general_table.py /usr/local/bin/build_general_table.py \
    && ln -s /opt/deletions_to_csv.py /usr/local/bin/deletions_to_csv.py \
    && ln -s /opt/profiler_parser.py /usr/local/bin/profiler_parser.py
