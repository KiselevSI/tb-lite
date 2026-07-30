#!/bin/bash
# build-docker-images.sh
#
# Сборка Docker-образов для запуска пайплайна в Docker-режиме.
# Кастомные образы собираются из Dockerfile-ов в containers/dockerfiles/
# Публичные образы скачиваются из Docker Hub.
#
# Использование:
#   cd /path/to/tb-lite
#   bash build-docker-images.sh
#
# После сборки запускайте пайплайн:
#   nextflow run main.nf -profile docker

set -euo pipefail

TARGET="${1:-all}"

# Контекст сборки — корень проекта (нужен для COPY assets/scripts/*)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

TEMP_DOCKER_CONFIG=""

cleanup() {
    if [[ -n "${TEMP_DOCKER_CONFIG:-}" && -d "$TEMP_DOCKER_CONFIG" ]]; then
        rm -rf "$TEMP_DOCKER_CONFIG"
    fi
}

configure_docker_auth() {
    local docker_config_dir docker_config_file creds_store helper

    docker_config_dir="${DOCKER_CONFIG:-$HOME/.docker}"
    docker_config_file="$docker_config_dir/config.json"

    if [[ ! -f "$docker_config_file" ]]; then
        return 0
    fi

    creds_store="$(sed -nE 's/^[[:space:]]*"credsStore"[[:space:]]*:[[:space:]]*"([^"]+)".*$/\1/p' "$docker_config_file" | head -n1)"
    if [[ -z "$creds_store" ]]; then
        return 0
    fi

    helper="docker-credential-$creds_store"
    if command -v "$helper" >/dev/null 2>&1; then
        return 0
    fi

    TEMP_DOCKER_CONFIG="$(mktemp -d /tmp/tb-lite-docker-config.XXXXXX)"
    printf '{\n  "auths": {}\n}\n' > "$TEMP_DOCKER_CONFIG/config.json"
    export DOCKER_CONFIG="$TEMP_DOCKER_CONFIG"

    echo "Предупреждение: helper $helper не найден; используется временный DOCKER_CONFIG без credsStore." >&2
}

trap cleanup EXIT
configure_docker_auth

echo "=== Сборка кастомных Docker-образов ==="
echo "Контекст сборки: $(pwd)"

build_spotyping_image() {
    docker build -t tb-lite/spotyping:2.0 -f containers/dockerfiles/spotyping.Dockerfile .
}

build_ann_table_image() {
    docker build -t tb-lite/ann-table:1.0 -f containers/dockerfiles/ann_table.Dockerfile .
}

build_report_images() {
    docker build -t tb-lite/build-table:1.1        -f containers/dockerfiles/build_table.Dockerfile .
    docker build -t tb-lite/tb-platform-tables:1.1 -f containers/dockerfiles/tb_platform_tables.Dockerfile .
}

build_all_images() {
    build_ann_table_image
    build_report_images
    docker build -t tb-lite/tb-mix:1.0          -f containers/dockerfiles/tb_mix.Dockerfile          .
    docker build -t tb-lite/rd-scanner:2.0      -f containers/dockerfiles/rd.Dockerfile              .
    build_spotyping_image
    docker build -t tb-lite/tblg:1.0            -f containers/dockerfiles/tblg.Dockerfile            .
}

case "$TARGET" in
    all)
        build_all_images
        ;;
    reports|reports-only|tables)
        build_report_images
        ;;
    ann-table|ann-table-only)
        build_ann_table_image
        ;;
    spotyping|spotyping-only)
        build_spotyping_image
        ;;
    *)
        echo "Неизвестная цель сборки: $TARGET" >&2
        echo "Использование: bash build-docker-images.sh [all|reports|ann-table|spotyping]" >&2
        exit 1
        ;;
esac

if [[ "$TARGET" == "all" ]]; then
    echo ""
    echo "=== Скачивание публичных Docker-образов ==="

    # QC / reporting
    docker pull staphb/fastqc:0.12.1
    docker pull multiqc/multiqc:v1.30

    # trimming / alignment / coverage
    docker pull community.wave.seqera.io/library/fastp:1.1.0--08aa7c5662a30d57
    docker pull community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9
    docker pull community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0
    docker pull community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5
    docker pull community.wave.seqera.io/library/mosdepth_htslib:0f58993cb6d93294

    # variant calling / bcftools / annotation
    docker pull quay.io/biocontainers/freebayes:1.3.10--hbefcdb2_0
    docker pull community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a
    docker pull community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6
    docker pull community.wave.seqera.io/library/snpeff:5.4.0c--e08ddc54579e82bd

    # genotyping / DR / IS6110
    docker pull staphb/tbprofiler:6.6.6
    docker pull quay.io/biocontainers/ismapper:2.0.2--pyhdfd78af_1

    # taxonomy (Kraken2 / Bracken)
    docker pull community.wave.seqera.io/library/kraken2_coreutils_pigz:920ecc6b96e2ba71
    docker pull community.wave.seqera.io/library/bracken:3.1--22a4e66ce04c5e01

    # SRA download
    docker pull quay.io/biocontainers/sra-tools:3.2.1--h4304569_1
    docker pull community.wave.seqera.io/library/sra-tools_pigz:4a694d823f6f7fcf

    # misc
    docker pull ubuntu:22.04
fi

echo ""
echo "=== Готово! ==="
echo ""
echo "Запуск пайплайна:"
echo "  nextflow run main.nf -profile docker"
