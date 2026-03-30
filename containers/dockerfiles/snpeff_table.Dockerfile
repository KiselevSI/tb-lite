# snpeff_table -- SnpEff + SnpSift + bcftools/samtools from local installation
# Build context: project root (tb-lite/)

FROM eclipse-temurin:21-jre-jammy

LABEL maintainer="TB-Lite pipeline"

ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    DEBIAN_FRONTEND=noninteractive

COPY snpEff_latest_core/snpEff /opt/snpEff

# Install bcftools, samtools and utilities needed by ANN_TABLE processes
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bcftools samtools \
        procps coreutils gawk grep sed findutils \
    && rm -rf /var/lib/apt/lists/* \
    # Make SnpSift available as a command
    && printf '#!/bin/bash\nexec java -jar /opt/snpEff/SnpSift.jar "$@"\n' > /usr/local/bin/snpSift \
    && chmod +x /usr/local/bin/snpSift \
    # Provide snpEff with the bundled local config and databases
    && printf '#!/bin/bash\nexec java -jar /opt/snpEff/snpEff.jar -c /opt/snpEff/snpEff.config "$@"\n' > /usr/local/bin/snpEff \
    && chmod +x /usr/local/bin/snpEff

ENV SNPEFF_HOME=/opt/snpEff
