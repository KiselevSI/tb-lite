#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Батчевый запуск TB-Lite pipeline
# Разбивает CSV на батчи, запускает каждый отдельно,
# после завершения очищает work/ для экономии места.
# ============================================================

# --- Значения по умолчанию ---
INPUT=""
PIPELINE_DIR=""
BATCH_SIZE=500
PROFILE="local"
OUTDIR=""
WORKDIR=""
RESUME_FROM=1
INPUT_MODE="auto"

usage() {
    cat <<EOF
Использование:
  $0 --input <samples.csv> [опции]

Опции:
  --input        CSV файл со всеми образцами (обязательный)
  --pipeline     Путь к папке с пайплайном TB-Lite (обязательный)
  --input-mode   Тип входа: auto, fastq, sra (по умолчанию: auto)
  --batch-size   Количество образцов в батче (по умолчанию: 500)
  --profile      Nextflow профиль (по умолчанию: local)
  --outdir       Папка для результатов (по умолчанию: <текущая_директория>/results)
  --workdir      Рабочая директория Nextflow (по умолчанию: <текущая_директория>/work)
  --resume-from  Начать с батча N (по умолчанию: 1)
  --help         Показать справку

Примеры:
  $0 --pipeline /home/zerg/git/tb-lite --input /data/all_50k.csv
  $0 --pipeline /home/zerg/git/tb-lite --input /data/all_50k.csv --batch-size 500 --resume-from 3

Особенности batch-режима:
  - запуск идёт с флагами --skip_final_reports --skip_kraken --skip_snp_matrix
  - MultiQC сохраняется отдельно для каждого батча в multiqc_batch_<N>/
EOF
    exit 0
}

# --- Разбор аргументов ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)       INPUT="$2";        shift 2 ;;
        --pipeline)    PIPELINE_DIR="$2"; shift 2 ;;
        --input-mode)  INPUT_MODE="$2";   shift 2 ;;
        --batch-size)  BATCH_SIZE="$2";  shift 2 ;;
        --profile)     PROFILE="$2";     shift 2 ;;
        --outdir)      OUTDIR="$2";      shift 2 ;;
        --workdir)     WORKDIR="$2";     shift 2 ;;
        --resume-from) RESUME_FROM="$2"; shift 2 ;;
        --help)        usage ;;
        *) echo "Неизвестный аргумент: $1"; usage ;;
    esac
done

if [[ -z "$INPUT" ]]; then
    echo "Ошибка: --input обязателен"
    usage
fi

if [[ -z "$PIPELINE_DIR" ]]; then
    echo "Ошибка: --pipeline обязателен"
    usage
fi

if [[ ! -f "$INPUT" ]]; then
    echo "Ошибка: файл $INPUT не найден"
    exit 1
fi

if [[ ! -f "${PIPELINE_DIR}/main.nf" ]]; then
    echo "Ошибка: main.nf не найден в ${PIPELINE_DIR}"
    exit 1
fi

case "$INPUT_MODE" in
    auto|fastq|sra) ;;
    *)
        echo "Ошибка: --input-mode должен быть одним из: auto, fastq, sra"
        exit 1
        ;;
esac

detect_input_mode() {
    local input_file="$1"
    local first_nonempty

    if [[ "$INPUT_MODE" != "auto" ]]; then
        echo "$INPUT_MODE"
        return 0
    fi

    first_nonempty="$(grep -m1 '[^[:space:]]' "$input_file" || true)"
    if [[ -z "$first_nonempty" ]]; then
        echo "Ошибка: входной файл $input_file пуст" >&2
        exit 1
    fi

    if [[ "$first_nonempty" == *","* ]]; then
        echo "fastq"
    else
        echo "sra"
    fi
}

# Преобразуем в абсолютные пути
INPUT="$(cd "$(dirname "$INPUT")" && pwd)/$(basename "$INPUT")"
PIPELINE_DIR="$(cd "$PIPELINE_DIR" && pwd)"
[[ -z "$OUTDIR" ]] && OUTDIR="$(pwd)/results"
[[ -z "$WORKDIR" ]] && WORKDIR="$(pwd)/work"
OUTDIR="$(mkdir -p "$OUTDIR" && cd "$OUTDIR" && pwd)"
LOG_FILE="$(pwd)/batches.log"
BATCH_DIR="$(pwd)/.batches"
INPUT_MODE="$(detect_input_mode "$INPUT")"
if [[ "$INPUT_MODE" == "fastq" ]]; then
    BATCH_EXT="csv"
    NF_INPUT_FLAG="--input"
else
    BATCH_EXT="txt"
    NF_INPUT_FLAG="--sra_ids"
fi

# --- Подготовка ---
HEADER=""
if [[ "$INPUT_MODE" == "fastq" ]]; then
    HEADER=$(head -1 "$INPUT")
    TOTAL_SAMPLES=$(tail -n +2 "$INPUT" | grep -c '[^[:space:]]' || true)
else
    TOTAL_SAMPLES=$(grep -c '[^[:space:]]' "$INPUT" || true)
fi
TOTAL_BATCHES=$(( (TOTAL_SAMPLES + BATCH_SIZE - 1) / BATCH_SIZE ))

echo "============================================"
echo "TB-Lite батчевый запуск"
echo "============================================"
echo "Входной файл:  $INPUT"
echo "Образцов:      $TOTAL_SAMPLES"
echo "Размер батча:  $BATCH_SIZE"
echo "Всего батчей:  $TOTAL_BATCHES"
echo "Начать с:      $RESUME_FROM"
echo "Режим входа:   $INPUT_MODE"
echo "Профиль:       $PROFILE"
echo "Результаты:    $OUTDIR"
echo "Work dir:      $WORKDIR"
echo "Пайплайн:      $PIPELINE_DIR"
echo "============================================"

mkdir -p "$BATCH_DIR"
mkdir -p "$OUTDIR"
mkdir -p "$WORKDIR"
rm -f "${BATCH_DIR}"/batch_*.csv "${BATCH_DIR}"/batch_*.txt

archive_batch_multiqc() {
    local batch_num="$1"
    local src_dir="${OUTDIR}/multiqc"
    local dst_dir="${OUTDIR}/multiqc_batch_${batch_num}"

    if [[ ! -d "$src_dir" ]]; then
        return 0
    fi

    if [[ -e "$dst_dir" ]]; then
        local backup_dir="${dst_dir}.prev_$(date '+%Y%m%d_%H%M%S')"
        echo "  Найден существующий ${dst_dir}, переношу его в ${backup_dir}"
        mv "$dst_dir" "$backup_dir"
    fi

    mv "$src_dir" "$dst_dir"
    echo "  MultiQC сохранён в ${dst_dir}"
}

# --- Разбиение входа на батчи ---
echo "Разбиваю входной файл на батчи..."

BATCH_NUM=0
LINE_NUM=0

if [[ "$INPUT_MODE" == "fastq" ]]; then
    tail -n +2 "$INPUT" | grep '[^[:space:]]' | while IFS= read -r line; do
        if (( LINE_NUM % BATCH_SIZE == 0 )); then
            BATCH_NUM=$(( LINE_NUM / BATCH_SIZE + 1 ))
            BATCH_FILE="${BATCH_DIR}/batch_${BATCH_NUM}.${BATCH_EXT}"
            echo "$HEADER" > "$BATCH_FILE"
        fi
        BATCH_FILE="${BATCH_DIR}/batch_$(( LINE_NUM / BATCH_SIZE + 1 )).${BATCH_EXT}"
        echo "$line" >> "$BATCH_FILE"
        LINE_NUM=$(( LINE_NUM + 1 ))
    done
else
    grep '[^[:space:]]' "$INPUT" | while IFS= read -r line; do
        if (( LINE_NUM % BATCH_SIZE == 0 )); then
            BATCH_NUM=$(( LINE_NUM / BATCH_SIZE + 1 ))
            BATCH_FILE="${BATCH_DIR}/batch_${BATCH_NUM}.${BATCH_EXT}"
            : > "$BATCH_FILE"
        fi
        BATCH_FILE="${BATCH_DIR}/batch_$(( LINE_NUM / BATCH_SIZE + 1 )).${BATCH_EXT}"
        echo "$line" >> "$BATCH_FILE"
        LINE_NUM=$(( LINE_NUM + 1 ))
    done
fi

echo "Создано $TOTAL_BATCHES батчей в $BATCH_DIR/"

# --- Запуск батчей ---
for (( i = RESUME_FROM; i <= TOTAL_BATCHES; i++ )); do
    BATCH_FILE="${BATCH_DIR}/batch_${i}.${BATCH_EXT}"
    if [[ "$INPUT_MODE" == "fastq" ]]; then
        BATCH_SAMPLES=$(( $(wc -l < "$BATCH_FILE") - 1 ))
    else
        BATCH_SAMPLES=$(wc -l < "$BATCH_FILE")
    fi

    echo ""
    echo "[batch ${i}/${TOTAL_BATCHES}] Запуск ${BATCH_SAMPLES} образцов..."

    # Определяем нужен ли -resume (для первого батча при --resume-from)
    NF_RESUME=""
    if (( i == RESUME_FROM )) && [[ -d "$WORKDIR" ]] && [[ -n "$(ls -A "$WORKDIR" 2>/dev/null)" ]]; then
        echo "  Обнаружен work/ — используется -resume"
        NF_RESUME="-resume"
    fi

    # Запуск Nextflow
    if nextflow run "${PIPELINE_DIR}/main.nf" \
        "$NF_INPUT_FLAG" "$BATCH_FILE" \
        --outdir "$OUTDIR" \
        --skip_final_reports \
        --skip_kraken \
        --skip_snp_matrix \
        -w "$WORKDIR" \
        -profile "$PROFILE" \
        $NF_RESUME; then

        echo "[batch ${i}/${TOTAL_BATCHES}] Завершён успешно"
        echo "batch_${i} OK $(date '+%Y-%m-%d %H:%M:%S') samples=${BATCH_SAMPLES}" >> "$LOG_FILE"

        archive_batch_multiqc "$i"

        # Очистка work/
        echo "  Очищаю work/..."
        rm -rf "${WORKDIR:?}"/*
        echo "  work/ очищен"
    else
        EXIT_CODE=$?
        echo ""
        echo "============================================"
        echo "ОШИБКА: batch ${i} завершился с кодом ${EXIT_CODE}"
        echo "============================================"
        echo "Для продолжения выполните:"
        echo "  $0 --input $INPUT --pipeline $PIPELINE_DIR --input-mode $INPUT_MODE --batch-size $BATCH_SIZE --profile $PROFILE --outdir $OUTDIR --workdir $WORKDIR --resume-from $i"
        echo ""
        echo "work/ сохранён для возможности -resume"
        exit "$EXIT_CODE"
    fi
done

echo ""
echo "============================================"
echo "Все $TOTAL_BATCHES батчей завершены!"
echo "Результаты в: $OUTDIR"
echo "Лог: $LOG_FILE"
echo "============================================"
