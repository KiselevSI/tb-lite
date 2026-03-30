# mosdepth -- Fast BAM/CRAM depth calculation
# Converted from: containers/def/mosdpt.def
# Build context: project root (tb-lite/)

FROM ubuntu:22.04

LABEL maintainer="TB-Lite pipeline"

ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PATH=/usr/bin:$PATH \
    DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       wget \
       ca-certificates locales \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    # Download mosdepth binary
    && cd /opt \
    && wget https://github.com/brentp/mosdepth/releases/download/v0.3.11/mosdepth \
    && chmod +x mosdepth \
    && ln -s /opt/mosdepth /usr/bin/mosdepth
