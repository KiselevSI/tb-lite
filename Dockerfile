###############################################################################
#  TB‑Lite pipeline — all‑in‑one Docker image (Java 17, Nextflow, Apptainer)
#  Build:   docker build -t tb-lite:2025-07 .
#  Run:     docker run --rm -it --privileged \
#                   -v $(pwd)/data:/workspace/data \
#                   -v $(pwd)/results:/workspace/results \
#                   tb-lite:2025-07 \
#                   run main.nf -with-singularity -resume
###############################################################################

############################ 1️⃣ Base image with Java ##########################
FROM eclipse-temurin:17-jre-jammy AS base

LABEL maintainer="you@example.org" \
      org.opencontainers.image.source="https://github.com/yourrepo/tb-lite"

ENV DEBIAN_FRONTEND=noninteractive \
    NXF_VER=25.04.6

############################ 2️⃣ System + build deps ###########################
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl git procps build-essential \
        squashfs-tools libseccomp-dev uuid-dev \
        fuse-overlayfs \
        tini \
    && rm -rf /var/lib/apt/lists/*            /* keep image slim */ :contentReference[oaicite:3]{index=3}

############################ 3️⃣ Nextflow #####################################
RUN curl -s https://get.nextflow.io | bash \
 && mv nextflow /usr/local/bin/ \
 && chmod +x /usr/local/bin/nextflow    

############################ 4️⃣ Apptainer #####################################
# 1.4.x поддерживает fuse-overlayfs без SUID; ставим бинарный .deb
RUN curl -L https://github.com/apptainer/apptainer/releases/download/v1.4.1/apptainer_1.4.1_amd64.deb \
      -o /tmp/apptainer.deb \
 && apt-get update && apt-get install -y /tmp/apptainer.deb \
 && rm /tmp/apptainer.deb 

############################ 5️⃣ Pipeline code #################################
WORKDIR /workspace

# Копируем всё, кроме файлов, исключённых .dockerignore  (FASTQ/BAM  и т.д.) :contentReference[oaicite:6]{index=6}
COPY main.nf nextflow.config scripts/ ref/ db/ IS6110_db/
COPY containers/ ./containers          

# Nextflow будет искать SIF‑файлы именно здесь
ENV NXF_SINGULARITY_CACHEDIR=/workspace/containers 

############################ 6️⃣ Entrypoint ####################################
# tini → корректный обработчик сигналов; Nextflow — основной процесс
ENTRYPOINT ["/usr/bin/tini","--","nextflow"]
CMD ["run","main.nf","-with-singularity","-resume"]
