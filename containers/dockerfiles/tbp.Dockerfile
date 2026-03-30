FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/conda/bin:$PATH"

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget bzip2 ca-certificates curl && \
    rm -rf /var/lib/apt/lists/*

RUN wget -qO /tmp/miniforge.sh \
    "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh" && \
    bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm /tmp/miniforge.sh

RUN conda create -y -n tbprofiler -c conda-forge -c bioconda \
    python=3.11 \
    "setuptools<82" \
    tb-profiler=6.6.6 \
    pathogen-profiler && \
    conda clean -afy

ENV PATH="/opt/conda/envs/tbprofiler/bin:$PATH"

COPY assets/h37rv.fa /ref/h37rv.fa

RUN tb-profiler version && \
    python -c "import pathogenprofiler; print(pathogenprofiler.__version__)"

# пока временно отключи
RUN tb-profiler update_tbdb --match_ref /ref/h37rv.fa