# freebayes -- Variant caller with bcftools and tabix
# Converted from: containers/def/freebayes.def
# Build context: project root (tb-lite/)

FROM ubuntu:24.04

LABEL maintainer="TB-Lite pipeline"

ENV DEBIAN_FRONTEND=noninteractive \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    PATH=/opt:$PATH \
    PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install -y --no-install-recommends \
      bcftools wget ca-certificates tabix \
    && rm -rf /var/lib/apt/lists/* \
    # Install freebayes static binary
    && mkdir -p /opt/freebayes \
    && cd /opt/freebayes \
    && wget https://github.com/freebayes/freebayes/releases/download/v1.3.10/freebayes-1.3.10-linux-amd64-static.gz \
    && gunzip freebayes-1.3.10-linux-amd64-static.gz \
    && chmod +x freebayes-1.3.10-linux-amd64-static \
    && mv freebayes-1.3.10-linux-amd64-static /usr/local/bin/freebayes
