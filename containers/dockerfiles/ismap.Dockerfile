# ismap -- IS_mapper for insertion sequence detection
# Converted from: containers/def/ismap.def
# Build context: project root (tb-lite/)

FROM python:3.8-slim

LABEL Author="Your Name" \
      Version="1.0" \
      maintainer="TB-Lite pipeline"

ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive

RUN set -e \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
       build-essential \
       git \
       bwa \
       samtools \
       bedtools \
       ncbi-blast+ procps \
    && rm -rf /var/lib/apt/lists/* \
    && pip install setuptools biopython pysam \
    && pip install git+https://github.com/jhawkey/IS_mapper
