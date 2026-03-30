# spotyping -- SpoTyping-v2.0 + BLAST via micromamba (Python 2.7)
# Converted from: containers/def/spotyping.def
# Build context: project root (tb-lite/)

FROM debian:bookworm-slim

LABEL Author="TB-Lite pipeline" \
      Version="2.2" \
      Description="SpoTyping-v2.0 + BLAST 2.12 (via micromamba)"

ENV LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8 \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/envs/spotenv/bin:/opt/spotyping:$PATH \
    DEBIAN_FRONTEND=noninteractive

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN set -eu \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
       locales ca-certificates wget curl bzip2 procps \
    && rm -rf /var/lib/apt/lists/*

RUN set -eu \
    && mkdir -p /opt/micromamba "$MAMBA_ROOT_PREFIX" \
    && curl -fL --retry 5 --retry-delay 2 --retry-all-errors \
       -o /tmp/micromamba.tar.bz2 \
       https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64.tar.bz2 \
    && bzip2 -t /tmp/micromamba.tar.bz2 \
    && tar -C /opt/micromamba -xvjf /tmp/micromamba.tar.bz2 \
       bin/micromamba --strip-components=1 \
    && rm -f /tmp/micromamba.tar.bz2 \
    && ln -s /opt/micromamba/micromamba /usr/local/bin/micromamba

RUN set -eu \
    && micromamba create -y -n spotenv \
       -c conda-forge -c bioconda \
       python=2.7 spotyping \
    && micromamba clean --all --yes
