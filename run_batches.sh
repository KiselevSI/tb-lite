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
BENCHMARK=0
ORIG_ARGS=("$@")

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
  --benchmark    Сохранить трассировку Nextflow и метрики производительности
                 в <outdir>/benchmark/ (время, CPU, память, образцов в час)
  --help         Показать справку

Примеры:
  $0 --pipeline /home/zerg/git/tb-lite --input /data/all_50k.csv
  $0 --pipeline /home/zerg/git/tb-lite --input /data/all_50k.csv --batch-size 500 --resume-from 3
  $0 --pipeline /home/zerg/git/tb-lite --input /data/all_50k.csv --with-kraken --kraken2_db /data/kraken_db

Особенности batch-режима:
  - каждый батч запускается с --skip_final_reports --skip_multiqc --skip_snp_matrix
  - после последнего батча автоматически собирается один общий Reports/
  - с --benchmark для каждого батча пишутся trace/report/timeline Nextflow,
    а в конце — сводка bin/benchmark_summary.py (summary.md, per_*.tsv)
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
        --benchmark)   BENCHMARK=1;      shift ;;
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

mkdir -p "$BATCH_DIR"
mkdir -p "$OUTDIR"
mkdir -p "$WORKDIR"
mkdir -p "${OUTDIR}/batch_reports/filter"
BENCH_DIR="${OUTDIR}/benchmark"

count_existing_batches() {
    local max_batch=0
    local batch_file
    local batch_name
    local batch_num

    shopt -s nullglob
    for batch_file in "${BATCH_DIR}"/batch_*.${BATCH_EXT}; do
        batch_name="${batch_file##*/}"
        batch_num="${batch_name#batch_}"
        batch_num="${batch_num%.${BATCH_EXT}}"
        if [[ "$batch_num" =~ ^[0-9]+$ ]] && (( batch_num > max_batch )); then
            max_batch="$batch_num"
        fi
    done
    shopt -u nullglob

    echo "$max_batch"
}

count_samples_in_batch_file() {
    local batch_file="$1"
    if [[ "$INPUT_MODE" == "fastq" ]]; then
        local lines
        lines=$(wc -l < "$batch_file")
        echo $(( lines > 0 ? lines - 1 : 0 ))
    else
        wc -l < "$batch_file"
    fi
}

count_samples_in_existing_batches() {
    local total=0
    local i
    local batch_file

    for (( i = 1; i <= TOTAL_BATCHES; i++ )); do
        batch_file="${BATCH_DIR}/batch_${i}.${BATCH_EXT}"
        if [[ ! -f "$batch_file" ]]; then
            echo "Ошибка: существующая разбивка неполная — нет ${batch_file}" >&2
            exit 1
        fi
        total=$(( total + $(count_samples_in_batch_file "$batch_file") ))
    done

    echo "$total"
}

create_batches_from_input() {
    echo "Создаю новую разбивку входного файла на батчи..."

    local header=""
    local batch_num=0
    local line_num=0
    local batch_file

    if [[ "$INPUT_MODE" == "fastq" ]]; then
        header=$(head -1 "$INPUT")
        tail -n +2 "$INPUT" | grep '[^[:space:]]' | while IFS= read -r line; do
            if (( line_num % BATCH_SIZE == 0 )); then
                batch_num=$(( line_num / BATCH_SIZE + 1 ))
                batch_file="${BATCH_DIR}/batch_${batch_num}.${BATCH_EXT}"
                echo "$header" > "$batch_file"
            fi
            batch_file="${BATCH_DIR}/batch_$(( line_num / BATCH_SIZE + 1 )).${BATCH_EXT}"
            echo "$line" >> "$batch_file"
            line_num=$(( line_num + 1 ))
        done
    else
        grep '[^[:space:]]' "$INPUT" | while IFS= read -r line; do
            if (( line_num % BATCH_SIZE == 0 )); then
                batch_num=$(( line_num / BATCH_SIZE + 1 ))
                batch_file="${BATCH_DIR}/batch_${batch_num}.${BATCH_EXT}"
                : > "$batch_file"
            fi
            batch_file="${BATCH_DIR}/batch_$(( line_num / BATCH_SIZE + 1 )).${BATCH_EXT}"
            echo "$line" >> "$batch_file"
            line_num=$(( line_num + 1 ))
        done
    fi
}

# --- Подготовка ---
TOTAL_BATCHES="$(count_existing_batches)"
if (( TOTAL_BATCHES > 0 )); then
    TOTAL_SAMPLES="$(count_samples_in_existing_batches)"
    echo "Использую существующую разбивку: $TOTAL_BATCHES батчей в $BATCH_DIR/"
else
    if [[ "$INPUT_MODE" == "fastq" ]]; then
        TOTAL_SAMPLES=$(tail -n +2 "$INPUT" | grep -c '[^[:space:]]' || true)
    else
        TOTAL_SAMPLES=$(grep -c '[^[:space:]]' "$INPUT" || true)
    fi
    TOTAL_BATCHES=$(( (TOTAL_SAMPLES + BATCH_SIZE - 1) / BATCH_SIZE ))
    create_batches_from_input
    TOTAL_BATCHES="$(count_existing_batches)"
    echo "Создано $TOTAL_BATCHES батчей в $BATCH_DIR/"
fi

if (( RESUME_FROM < 1 || RESUME_FROM > TOTAL_BATCHES )); then
    echo "Ошибка: --resume-from должен быть в диапазоне 1..${TOTAL_BATCHES}"
    exit 1
fi

if [[ ! -f "${BATCH_DIR}/batch_${RESUME_FROM}.${BATCH_EXT}" ]]; then
    echo "Ошибка: не найден ${BATCH_DIR}/batch_${RESUME_FROM}.${BATCH_EXT}"
    exit 1
fi

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
if (( BENCHMARK )); then
    echo "Benchmark:     $BENCH_DIR"
fi
echo "============================================"

# --- Бенчмарк ---
# Одна строка "ключ<TAB>значение"; табы и переводы строк в значении схлопываются.
env_row() {
    local key="$1"
    shift
    local value
    value="$("$@" 2>/dev/null | tr '\t\n' '  ' | sed 's/  */ /g; s/^ //; s/ $//' || true)"
    printf '%s\t%s\n' "$key" "${value:-NA}"
}

write_benchmark_environment() {
    local env_file="${BENCH_DIR}/environment.tsv"
    local engine_cmd=(true)

    case "$PROFILE" in
        docker)      engine_cmd=(docker --version) ;;
        singularity) engine_cmd=(singularity --version) ;;
        conda)       engine_cmd=(conda --version) ;;
    esac

    {
        printf 'key\tvalue\n'
        env_row date              date '+%Y-%m-%d %H:%M:%S %z'
        env_row hostname          hostname
        env_row os                sh -c '. /etc/os-release && echo "$PRETTY_NAME"'
        env_row kernel            uname -r
        env_row cpu_model         sh -c "lscpu | sed -n 's/^Model name:[[:space:]]*//p' | head -1"
        env_row cpu_sockets       sh -c "lscpu | sed -n 's/^Socket(s):[[:space:]]*//p'"
        env_row cpu_logical       nproc
        env_row mem_total_gb      env LC_ALL=C awk '/^MemTotal:/ { printf "%.1f", $2 / 1024 / 1024 }' /proc/meminfo
        env_row workdir_fs        sh -c "df -hT '$WORKDIR' | tail -1"
        env_row outdir_fs         sh -c "df -hT '$OUTDIR' | tail -1"
        env_row block_devices     sh -c "lsblk -d -n -o NAME,ROTA,SIZE,MODEL | sed 's/\$/;/'"
        env_row nextflow_version  sh -c "nextflow -version | sed -n 's/.*version \\([0-9][^ ]*\\).*/\\1/p'"
        env_row container_engine  "${engine_cmd[@]}"
        env_row pipeline_commit   git -C "$PIPELINE_DIR" rev-parse HEAD
        env_row pipeline_dirty    sh -c "git -C '$PIPELINE_DIR' status --porcelain --untracked-files=no | wc -l"
        env_row executor_cpus     sh -c "nextflow config -flat '$PIPELINE_DIR' | sed -n \"s/^executor.cpus = //p\""
        env_row profile           echo "$PROFILE"
        env_row input_mode        echo "$INPUT_MODE"
        env_row batch_size        echo "$BATCH_SIZE"
        env_row total_samples     echo "$TOTAL_SAMPLES"
        env_row total_batches     echo "$TOTAL_BATCHES"
        env_row kraken            echo "$WITH_KRAKEN"
        env_row command           echo "$0 ${ORIG_ARGS[*]}"
    } > "$env_file"
    echo "  Окружение для бенчмарка: ${env_file}"
}

# Добавляет в массив NF_CMD (по имени) флаги трассировки для запуска <label>.
# Трассировка прошлой попытки того же запуска (после падения) переносится в
# attempts/, чтобы benchmark_summary.py мог взять метрики задач, которые при
# -resume придут как CACHED.
add_benchmark_flags() {
    local -n cmd_ref="$1"
    local label="$2"
    local run_dir="${BENCH_DIR}/${label}"

    if [[ -f "${run_dir}/trace.tsv" ]]; then
        mkdir -p "${BENCH_DIR}/attempts"
        mv "$run_dir" "${BENCH_DIR}/attempts/${label}.$(date +%s)"
    fi
    mkdir -p "$run_dir"

    cmd_ref+=(
        -c "${PIPELINE_DIR}/conf/benchmark.config"
        -with-trace "${run_dir}/trace.tsv"
        -with-report "${run_dir}/report.html"
        -with-timeline "${run_dir}/timeline.html"
    )
}

# Строка в batches.tsv: label, samples, start, end, wall_sec, exit_code, work_bytes
record_benchmark_run() {
    local label="$1" samples="$2" start="$3" end="$4" exit_code="$5" work_dir="$6"
    local runs_file="${BENCH_DIR}/batches.tsv"
    local work_bytes="NA"

    if [[ -d "$work_dir" ]]; then
        work_bytes="$(du -sb "$work_dir" 2>/dev/null | cut -f1 || true)"
        work_bytes="${work_bytes:-NA}"
    fi
    if [[ ! -f "$runs_file" ]]; then
        printf 'run\tsamples\tstart_epoch\tend_epoch\twall_sec\texit_code\twork_bytes\n' > "$runs_file"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$label" "$samples" "$start" "$end" "$(( end - start ))" "$exit_code" "$work_bytes" >> "$runs_file"
}

run_benchmark_summary() {
    echo ""
    echo "Считаю сводку бенчмарка..."
    if python3 "${PIPELINE_DIR}/bin/benchmark_summary.py" --bench-dir "$BENCH_DIR" --outdir "$OUTDIR"; then
        echo "  Сводка: ${BENCH_DIR}/summary.md"
    else
        echo "  ПРЕДУПРЕЖДЕНИЕ: сводка бенчмарка не построена; запустите вручную:"
        echo "  python3 ${PIPELINE_DIR}/bin/benchmark_summary.py --bench-dir $BENCH_DIR --outdir $OUTDIR"
    fi
}

if (( BENCHMARK )); then
    mkdir -p "$BENCH_DIR"
    write_benchmark_environment
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

merge_invalid_fastqs() {
    local output_dir="${OUTDIR}/Reports/general"
    local output_file="${output_dir}/bad_reads_invalid_fastq.txt"
    local files=( "${OUTDIR}/batch_reports/filter"/bad_reads_invalid_fastq.batch_*.txt )

    mkdir -p "$output_dir"

    if [[ ! -e "${files[0]}" ]]; then
        return 0
    fi

    awk 'FNR == 1 && ++seen > 1 { next } { print }' "${files[@]}" > "$output_file"
    echo "  Сводный bad_reads_invalid_fastq сохранён в ${output_file}"
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
        --skip_multiqc
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

    if (( BENCHMARK )); then add_benchmark_flags nf_cmd final_reports; fi

    echo ""
    echo "Собираю общий Reports/..."
    local start_ts end_ts exit_code=0
    start_ts=$(date +%s)
    "${nf_cmd[@]}" || exit_code=$?
    end_ts=$(date +%s)
    if (( BENCHMARK )); then
        record_benchmark_run final_reports "$TOTAL_SAMPLES" "$start_ts" "$end_ts" "$exit_code" "$reports_workdir"
    fi
    return "$exit_code"
}

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
    if (( BENCHMARK )); then add_benchmark_flags NF_CMD "batch_${i}"; fi

    BATCH_START=$(date +%s)
    NF_EXIT=0
    "${NF_CMD[@]}" || NF_EXIT=$?
    if (( BENCHMARK )); then
        # du до очистки work/ — это пиковый объём work/ для батча
        record_benchmark_run "batch_${i}" "$BATCH_SAMPLES" "$BATCH_START" "$(date +%s)" "$NF_EXIT" "$WORKDIR"
    fi

    if (( NF_EXIT == 0 )); then

        echo "[batch ${i}/${TOTAL_BATCHES}] Завершён успешно"
        echo "batch_${i} OK $(date '+%Y-%m-%d %H:%M:%S') samples=${BATCH_SAMPLES}" >> "$LOG_FILE"

        # Очистка work/
        echo "  Очищаю work/..."
        rm -rf "${WORKDIR:?}"/*
        echo "  work/ очищен"
    else
        EXIT_CODE=$NF_EXIT
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
        if (( BENCHMARK )); then echo -n " --benchmark"; fi
        echo ""
        echo ""
        echo "work/ сохранён для возможности -resume"
        exit "$EXIT_CODE"
    fi
done

merge_bad_reads
merge_invalid_fastqs
merge_unsupported_layouts
run_final_reports
if (( BENCHMARK )); then run_benchmark_summary; fi

echo ""
echo "============================================"
echo "Все $TOTAL_BATCHES батчей завершены!"
echo "Результаты в: $OUTDIR"
echo "Итоговые отчёты: ${OUTDIR}/Reports"
echo "Лог: $LOG_FILE"
if (( BENCHMARK )); then echo "Бенчмарк: ${BENCH_DIR}/summary.md"; fi
echo "============================================"
