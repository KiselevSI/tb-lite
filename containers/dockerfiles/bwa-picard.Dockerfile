# bwa-picard -- BWA aligner + Picard tools on Ubuntu 20.04
# Converted from: containers/def/bwa-picard.def
# Build context: project root (tb-lite/)

FROM ubuntu:20.04

LABEL maintainer="TB-Lite pipeline"

ENV DEBIAN_FRONTEND=noninteractive \
    PATH=/usr/local/bin:$PATH

# Pre-configure timezone to avoid interactive prompts
RUN echo 'tzdata tzdata/Areas select Etc'    | debconf-set-selections \
    && echo 'tzdata tzdata/Zones/Etc select UTC' | debconf-set-selections \
    && apt-get update && apt-get install -y \
       bwa \
       samtools \
       openjdk-17-jre-headless \
       curl \
    && apt-get clean && rm -rf /var/lib/apt/lists/* \
    # Download Picard.jar
    && curl -L -o /usr/local/bin/picard.jar \
       https://github.com/broadinstitute/picard/releases/download/3.4.0/picard.jar \
    && chmod +x /usr/local/bin/picard.jar

CMD ["java", "-jar", "/usr/local/bin/picard.jar"]
