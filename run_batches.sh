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
PROFILE="conda"
OUTDIR=""
WORKDIR=""
RESUME_FROM=1
INPUT_MODE="auto"
WITH_KRAKEN=0
KRAKEN2_DB=""
KRAKEN2_DB_LABEL=""
KRAKEN2_DB_2=""
KRAKEN2_DB_LABEL_2=""

usage() {
    cat <<EOF
Использование:
  $0 --input <samples.csv> [опции]

Опции:
  --input        CSV файл со всеми образцами (обязательный)
  --pipeline     Путь к папке с пайплайном TB-Lite (обязательный)
  --input-mode   Тип входа: auto, fastq, sra (по умолчанию: auto)
  --batch-size   Количество образцов в батче (по умолчанию: 500)
  --profile      Nextflow профиль: local, docker, singularity, conda (по умолчанию: conda)
  --outdir       Папка для результатов (по умолчанию: <текущая_директория>/results)
  --workdir      Рабочая директория Nextflow (по умолчанию: <текущая_директория>/work)
  --resume-from  Начать с батча N (по умолчанию: 1)
  --with-kraken  Включить Kraken/Bracken в батчах
  --kraken2_db   Путь к первой Kraken2 БД
  --kraken2_db_label    Лейбл первой Kraken2 БД
  --kraken2_db_2        Путь ко второй Kraken2 БД
  --kraken2_db_label_2  Лейбл второй Kraken2 БД
  --help         Показать справку

Примеры:
  $0 --pipeline /home/zerg/git/tb-lite --input /data/all_50k.csv
  $0 --pipeline /home/zerg/git/tb-lite --input /data/all_50k.csv --batch-size 500 --resume-from 3
  $0 --pipeline /home/zerg/git/tb-lite --input /data/all_50k.csv --with-kraken --kraken2_db /data/kraken_db

Особенности batch-режима:
  - каждый батч запускается с --skip_final_reports --skip_multiqc --skip_snp_matrix
  - после последнего батча автоматически собираются один общий Reports/ и один общий multiqc/
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
        --with-kraken) WITH_KRAKEN=1;    shift ;;
        --kraken2_db)  KRAKEN2_DB="$2";  shift 2 ;;
        --kraken2_db_label) KRAKEN2_DB_LABEL="$2"; shift 2 ;;
        --kraken2_db_2) KRAKEN2_DB_2="$2"; shift 2 ;;
        --kraken2_db_label_2) KRAKEN2_DB_LABEL_2="$2"; shift 2 ;;
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

case "$PROFILE" in
    local|docker|singularity|conda) ;;
    *)
        echo "Ошибка: --profile должен быть одним из: local, docker, singularity, conda"
        exit 1
        ;;
esac

if [[ -n "$KRAKEN2_DB" || -n "$KRAKEN2_DB_LABEL" || -n "$KRAKEN2_DB_2" || -n "$KRAKEN2_DB_LABEL_2" ]]; then
    WITH_KRAKEN=1
fi

if (( WITH_KRAKEN )) && [[ -z "$KRAKEN2_DB" ]]; then
    echo "Ошибка: для Kraken укажите --kraken2_db"
    exit 1
fi

if [[ -n "$KRAKEN2_DB" && ! -d "$KRAKEN2_DB" ]]; then
    echo "Ошибка: Kraken2 БД не найдена: $KRAKEN2_DB"
    exit 1
fi

if [[ -n "$KRAKEN2_DB_2" && ! -d "$KRAKEN2_DB_2" ]]; then
    echo "Ошибка: вторая Kraken2 БД не найдена: $KRAKEN2_DB_2"
    exit 1
fi

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
if [[ -n "$KRAKEN2_DB" ]]; then
    KRAKEN2_DB="$(cd "$(dirname "$KRAKEN2_DB")" && pwd)/$(basename "$KRAKEN2_DB")"
fi
if [[ -n "$KRAKEN2_DB_2" ]]; then
    KRAKEN2_DB_2="$(cd "$(dirname "$KRAKEN2_DB_2")" && pwd)/$(basename "$KRAKEN2_DB_2")"
fi
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
if (( WITH_KRAKEN )); then
    echo "Kraken:        enabled"
    echo "  DB1:         $KRAKEN2_DB"
    [[ -n "$KRAKEN2_DB_LABEL" ]] && echo "  DB1 label:   $KRAKEN2_DB_LABEL"
    [[ -n "$KRAKEN2_DB_2" ]] && echo "  DB2:         $KRAKEN2_DB_2"
    [[ -n "$KRAKEN2_DB_LABEL_2" ]] && echo "  DB2 label:   $KRAKEN2_DB_LABEL_2"
else
    echo "Kraken:        disabled"
fi
echo "============================================"

mkdir -p "$BATCH_DIR"
mkdir -p "$OUTDIR"
mkdir -p "$WORKDIR"
mkdir -p "${OUTDIR}/batch_reports/filter"
rm -f "${BATCH_DIR}"/batch_*.csv "${BATCH_DIR}"/batch_*.txt
if (( RESUME_FROM == 1 )); then
    rm -f "${OUTDIR}/batch_reports/filter"/bad_reads_low_coverage.batch_*.txt
    rm -f "${OUTDIR}/batch_reports/filter"/unsupported_sra_layout.batch_*.txt
fi

merge_bad_reads() {
    local output_dir="${OUTDIR}/Reports/general"
    local output_file="${output_dir}/bad_reads_low_coverage.txt"
    local files=( "${OUTDIR}/batch_reports/filter"/bad_reads_low_coverage.batch_*.txt )

    mkdir -p "$output_dir"

    if [[ ! -e "${files[0]}" ]]; then
        return 0
    fi

    awk 'FNR == 1 && ++seen > 1 { next } { print }' "${files[@]}" > "$output_file"
    echo "  Сводный bad_reads сохранён в ${output_file}"
}

merge_unsupported_layouts() {
    local output_dir="${OUTDIR}/Reports/general"
    local output_file="${output_dir}/unsupported_sra_layout.txt"
    local files=( "${OUTDIR}/batch_reports/filter"/unsupported_sra_layout.batch_*.txt )

    mkdir -p "$output_dir"

    if [[ ! -e "${files[0]}" ]]; then
        return 0
    fi

    awk 'FNR == 1 && ++seen > 1 { next } { print }' "${files[@]}" > "$output_file"
    echo "  Сводный unsupported_sra_layout сохранён в ${output_file}"
}

run_final_reports() {
    local reports_workdir="${WORKDIR}/_batch_reports"
    local nf_cmd=(
        nextflow run "${PIPELINE_DIR}/batch_reports.nf"
        --outdir "$OUTDIR"
        -w "$reports_workdir"
        -resume
    )

    if [[ "$PROFILE" != "local" ]]; then
        nf_cmd+=(-profile "$PROFILE")
    fi

    if (( WITH_KRAKEN )); then
        nf_cmd+=(--kraken2_db "$KRAKEN2_DB")
        [[ -n "$KRAKEN2_DB_LABEL" ]] && nf_cmd+=(--kraken2_db_label "$KRAKEN2_DB_LABEL")
        [[ -n "$KRAKEN2_DB_2" ]] && nf_cmd+=(--kraken2_db_2 "$KRAKEN2_DB_2")
        [[ -n "$KRAKEN2_DB_LABEL_2" ]] && nf_cmd+=(--kraken2_db_label_2 "$KRAKEN2_DB_LABEL_2")
    else
        nf_cmd+=(--skip_kraken)
    fi

    echo ""
    echo "Собираю общий Reports/ и multiqc/..."
    "${nf_cmd[@]}"
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

    NF_CMD=(
        nextflow run "${PIPELINE_DIR}/main.nf"
        "$NF_INPUT_FLAG" "$BATCH_FILE"
        --outdir "$OUTDIR"
        --batch_tag "batch_${i}"
        --skip_final_reports
        --skip_multiqc
        --skip_snp_matrix
        -w "$WORKDIR"
    )

    if [[ "$PROFILE" != "local" ]]; then
        NF_CMD+=(-profile "$PROFILE")
    fi

    if (( WITH_KRAKEN )); then
        NF_CMD+=(--kraken2_db "$KRAKEN2_DB")
        [[ -n "$KRAKEN2_DB_LABEL" ]] && NF_CMD+=(--kraken2_db_label "$KRAKEN2_DB_LABEL")
        [[ -n "$KRAKEN2_DB_2" ]] && NF_CMD+=(--kraken2_db_2 "$KRAKEN2_DB_2")
        [[ -n "$KRAKEN2_DB_LABEL_2" ]] && NF_CMD+=(--kraken2_db_label_2 "$KRAKEN2_DB_LABEL_2")
    else
        NF_CMD+=(--skip_kraken)
    fi

    [[ -n "$NF_RESUME" ]] && NF_CMD+=("$NF_RESUME")

    if "${NF_CMD[@]}"; then

        echo "[batch ${i}/${TOTAL_BATCHES}] Завершён успешно"
        echo "batch_${i} OK $(date '+%Y-%m-%d %H:%M:%S') samples=${BATCH_SAMPLES}" >> "$LOG_FILE"

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
        echo -n "  $0 --input $INPUT --pipeline $PIPELINE_DIR --input-mode $INPUT_MODE --batch-size $BATCH_SIZE --profile $PROFILE --outdir $OUTDIR --workdir $WORKDIR --resume-from $i"
        if (( WITH_KRAKEN )); then
            echo -n " --with-kraken --kraken2_db $KRAKEN2_DB"
            [[ -n "$KRAKEN2_DB_LABEL" ]] && echo -n " --kraken2_db_label $KRAKEN2_DB_LABEL"
            [[ -n "$KRAKEN2_DB_2" ]] && echo -n " --kraken2_db_2 $KRAKEN2_DB_2"
            [[ -n "$KRAKEN2_DB_LABEL_2" ]] && echo -n " --kraken2_db_label_2 $KRAKEN2_DB_LABEL_2"
        fi
        echo ""
        echo ""
        echo "work/ сохранён для возможности -resume"
        exit "$EXIT_CODE"
    fi
done

merge_bad_reads
merge_unsupported_layouts
run_final_reports

echo ""
echo "============================================"
echo "Все $TOTAL_BATCHES батчей завершены!"
echo "Результаты в: $OUTDIR"
echo "Итоговые отчёты: ${OUTDIR}/Reports"
echo "Итоговый MultiQC: ${OUTDIR}/multiqc"
echo "Лог: $LOG_FILE"
echo "============================================"
