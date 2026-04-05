# TB-Lite: Nextflow-пайплайн для геномного анализа M. tuberculosis

## Обзор

TB-Lite — это Nextflow-пайплайн для полного геномного анализа *Mycobacterium tuberculosis* (палочка Коха). Пайплайн принимает либо локальные gzipped FASTQ-файлы секвенирования, либо список SRA accession ID, и выполняет полный цикл анализа: от контроля качества сырых ридов до определения филогенетической линии, споликотипа, регионов делеций (RD) и предсказания лекарственной устойчивости.

**Референсный геном:** H37Rv (*M. tuberculosis*)

### Ключевые особенности

- **nf-core-совместимый каркас** — стандартные шаги QC/reporting и Kraken2/Bracken подключены через nf-core modules
- **Автоматическая фильтрация** — образцы с низким покрытием (< 10X медиана) автоматически исключаются из дальнейшего анализа
- **Детекция смешанных инфекций** — модуль TB-Mix выявляет образцы, содержащие несколько штаммов MTB
- **Мультиплатформенность** — стандартные runtime-профили `docker`, `singularity`, `conda`; по умолчанию включён Docker
- **Два режима входа** — `--input` для FASTQ samplesheet и `--sra_ids` для списка SRA accession
- **Опциональный Kraken2/Bracken** — отдельная ветка таксономической классификации с собственным флагом запуска
- **Когортный анализ** — при подаче > 1 образца автоматически строится объединённая таблица аннотированных вариантов
- **Поддержка paired-end и single-end** ридов (ISMapper работает только с парными)

---

## Структура проекта

```
tb-lite/
├── main.nf                    # Главный workflow (точка входа)
├── nextflow.config            # Конфигурация: параметры, профили, контейнеры
├── Dockerfile                 # Docker-образ для деплоя на сервере
├── build-containers.sh        # Скрипт сборки Singularity-контейнеров
├── run.csv                    # Пример FASTQ samplesheet
│
├── subworkflows/              # 8 подпроцессов (логические блоки)
│   ├── trimming.nf            # Тримминг ридов
│   ├── qc.nf                  # Контроль качества
│   ├── mapping.nf             # Картирование на референс
│   ├── filter.nf              # Фильтрация по качеству/покрытию
│   ├── call_varints.nf        # Вызов вариантов
│   ├── genotyping.nf          # Генотипирование и типирование
│   ├── reports.nf             # Отчёты
│   └── ann_table.nf           # Когортная таблица аннотаций
│
├── modules/
│   ├── nf-core/               # Готовые nf-core modules: fastp, fastqc, bwa, picard, samtools, freebayes, bcftools, multiqc, kraken2, bracken
│   └── local/                 # TB-специфичные и табличные модули: tb_mix, spotyping, rd, tblg, reports, ann_table post-processing
│
├── assets/                    # Референсные данные
│   ├── h37rv.fa               # Референсный геном H37Rv
│   ├── h37rv.gbk              # GenBank-файл H37Rv
│   ├── rd/RD.bed              # Регионы различий (RD)
│   ├── ismap/is6110.fasta     # Последовательность IS6110
│   ├── chr_name/chr.txt       # Маппинг имён хромосом
│   ├── tbmix/levels.tsv       # Уровни для детекции смешанных инфекций
│   ├── multiqc/multiqc_config.yaml
│   └── SNPEFF_ANNOTATION/h37rv_feature_table.txt
│
├── containers/                # Singularity-образы (.sif)
│   └── def/                   # Singularity definition files
│
├── data/                      # FASTQ-файлы (входные данные)
└── results/                   # Результаты анализа
```

---

## Входные данные

### 1. FASTQ samplesheet (`--input`)

CSV-файл с колонками:

| Колонка    | Описание |
|------------|----------|
| `sample`   | Имя образца |
| `fastq_1`  | Полный путь к `*.fastq.gz` или `*.fq.gz` для R1 или single-end |
| `fastq_2`  | Полный путь к `*.fastq.gz` или `*.fq.gz` для R2, пусто для single-end |

Совместимость со старым форматом `Sample,R1,R2,Layout` сохранена.

Важно: для `--input` поддерживаются только gzipped FASTQ (`*.fastq.gz` или `*.fq.gz`). Это тот же контракт, что и в `asm-lite`. Если у вас локальные файлы в виде `*.fastq`, их нужно предварительно сжать.

**Пример:**
```csv
sample,fastq_1,fastq_2
ERR123,/data/ERR123_1.fastq.gz,/data/ERR123_2.fastq.gz
SRR456,/data/SRR456.fastq.gz,
```

### 2. Список SRA accession (`--sra_ids`)

Текстовый файл: один accession на строку.

**Пример:**
```text
SRR32010433
ERR15166664
```

Пайплайн скачивает данные через nf-core modules `sratools/prefetch` и `sratools/fasterqdump`, автоматически определяет single-end/paired-end layout и публикует полученные FASTQ в `outdir/sra_fastq/`.

---

## Конвейер обработки

### Общая схема (8 этапов)

```
run.csv
    │
    ▼
┌─────────────────────────────────────┐
│ 1. TRIMMING (fastp)                 │  Тримминг адаптеров, polyG, качества
│    nf-core: fastp                  │
└─────────────┬───────────────────────┘
              │
    ┌─────────┴─────────┐
    ▼                   ▼
┌──────────┐    ┌──────────────────────┐
│ 2. QC    │    │ 3. MAPPING           │
│ (FastQC) │    │ nf-core bwa/mem      │
│          │    │ → picard/markduplicates │
└──────────┘    └─────────┬────────────┘
                          │
                          ▼
              ┌──────────────────────────────┐
              │ 4. FILTER                     │
              │  • TB-Mix (смешанная инфекция)│
              │  • Picard WGS/Align metrics   │
              │  • Samtools stats/flagstat     │
              │  Порог: median_cov >= 10X     │
              └─────────┬────────────────────┘
                        │ (только good samples)
                        ▼
              ┌──────────────────────────────┐
              │ 5. VARIANT CALLING            │
              │  Freebayes (ploidy=1)         │
              │  → BCFtools norm/filter       │
              │  → SnpEff аннотация           │
              │  → BCFtools stats             │
              └─────────┬────────────────────┘
                        │
              ┌─────────┴─────────┐
              ▼                   ▼
┌──────────────────┐  ┌──────────────────────┐
│ 6. GENOTYPING    │  │ 8. ANN_TABLE         │
│  (параллельно)   │  │  Когортный VCF       │
│  • SpoTyping     │  │  → SnpEff            │
│  • ISMapper      │  │  → SnpSift → TSV     │
│  • Mosdepth      │  └──────────────────────┘
│  • RD-анализ     │
│  • TBLG          │
│  • TB-Profiler   │
└────────┬─────────┘
         │
         ▼
┌──────────────────────────────┐
│ 7. REPORTS                    │
│  • MultiQC (HTML-дашборд)     │
│  • FINAL_TABLE.xlsx           │
└──────────────────────────────┘
```

### Как работает пайплайн: пошаговое описание

Пайплайн построен по принципу **последовательных этапов с точками ветвления**. Каждый этап получает данные из предыдущего через Nextflow-каналы, что обеспечивает параллельную обработку образцов.

**Общая логика потока данных:**

1. На вход подаётся CSV-файл с путями к FASTQ-файлам
2. Риды проходят тримминг (fastp) и далее **разделяются на два параллельных потока**: QC (FastQC) и картирование (BWA)
3. После картирования BAM-файлы проходят **этап фильтрации** — это ключевая точка пайплайна, где образцы делятся на «хорошие» и «плохие» по медианному покрытию (порог 10X)
4. Только «хорошие» образцы продолжают путь через вызов вариантов (Freebayes) и генотипирование
5. Генотипирование запускает **6 параллельных инструментов** для определения линии, споликотипа, делеций, покрытия и лекарственной устойчивости
6. Все результаты собираются в итоговые отчёты (MultiQC + FINAL_TABLE.xlsx)
7. При наличии > 1 образца и если не задан `--skip_snp_matrix`, запускается когортный анализ (ANN_TABLE)

**Ключевые точки ветвления:**

- **После TRIMMING** → данные идут в QC (параллельно) и MAPPING
- **После FILTER** → «хорошие» образцы передаются в CALLVAR и GENOTYPE; «плохие» записываются в `bad_reads_low_coverage.txt`
- **После CALLVAR** → VCF используется одновременно для GENOTYPE (TBLG-классификация), REPORTS (bcftools_stats) и, при необходимости, ANN_TABLE (когортное объединение)
- **GENOTYPE** → все 6 модулей запускаются **параллельно** друг другу

### Подробное описание этапов

#### 1. TRIMMING — Тримминг ридов

**Модули:** nf-core `fastp`

- **fastp** — обрезка адаптеров (автодетекция), polyG-хвостов, фильтрация по качеству
**Выход:** Триммированные FASTQ-файлы и JSON/HTML-отчёты fastp

#### 2. QC — Контроль качества

**Модуль:** nf-core `fastqc`

- Генерация отчётов FastQC для каждого образца
- Оценка качества секвенирования, GC-состава, длин ридов

#### 3. MAPPING — Картирование

**Модули:** nf-core `bwa/index`, `bwa/mem`, `samtools/faidx`, `picard/markduplicates`

- **BWA-mem** — выравнивание ридов на референс H37Rv
- **Picard MarkDuplicates** — маркировка PCR-дупликатов и построение BAM index

**Выход:** Отсортированные, дедуплицированные BAM-файлы

#### 4. FILTER — Фильтрация образцов

**Модули:** local `tb_mix`, nf-core `picard/collectwgsmetrics`, `picard/collectalignmentsummarymetrics`, `samtools/stats`, `samtools/flagstat`

Это **критический этап**, который разделяет образцы на два потока. Три модуля запускаются параллельно:

- **TB-Mix** (tb_mix.py) — детекция смешанных инфекций (несколько штаммов в образце). Использует BAM-файл и таблицу уровней (`levels.tsv`) для оценки гетерогенности аллельных частот
- **Picard CollectWgsMetrics** — рассчитывает медианное и среднее покрытие генома. Извлекает `mean_coverage` и `median_coverage` из Picard-отчёта
- **Picard CollectAlignmentSummaryMetrics** — метрики выравнивания
- **Samtools stats/flagstat** — базовая статистика BAM-файлов (общее число ридов, % картированных, % дупликатов)

**Логика фильтрации:**

Метрики каждого образца собираются в единый канал через `join()`. Затем фильтр отсеивает образцы по **медианному покрытию** (единственный активный фильтр):

- Медианное покрытие >= 30X (`min_median`) — **активный фильтр**
- % выравненных > 90% (`min_align_pct`) — **активный фильтр**

Образцы, не прошедшие фильтр → `bad_reads_low_coverage.txt` (с указанием причины). Прошедшие образцы передаются как `bam_good` и `trimmed_good` в последующие этапы

#### 5. VARIANT CALLING — Вызов вариантов

**Модули:** nf-core `freebayes`, `bcftools/index`, `bcftools/norm`, `bcftools/view`, `bcftools/stats`, `bcftools/annotate`, `snpeff/snpeff`

**Freebayes** (основной вариант-каллер):
- Гаплоидный режим (`--ploidy 1`) — TB не имеет диплоидного генома
- Минимальное качество картирования: 20
- Минимальное качество баз: 20
- Минимальное покрытие: 10X
- Минимальная частота альтернативного аллеля: 80%

**Постобработка:**
- BCFtools index + norm — индексация и нормализация VCF
- BCFtools view — фильтрация по `QUAL >= 10` и `INFO/DP >= 10`
- BCFtools annotate `--rename-chrs` — переименование хромосом по `chr.txt`
- SnpEff — аннотация вариантов

#### 6. GENOTYPING — Генотипирование

**Модули** (запускаются параллельно):

| Модуль | Инструмент | Что делает |
|--------|-----------|------------|
| `spotyping` | SpoTyping | Определяет споликотип (24-спейсерный бинарный код) |
| `ismapper` | nf-core `ismapper` | Картирование IS6110 (только для парных ридов) |
| `mosdepth` | nf-core `mosdepth` | Быстрый расчёт глубины покрытия |
| `rd` | rd.py | Анализ Region of Difference (RD) для линии |
| `tblg` | TBLG | Классификация линии (L1–L7) по VCF |
| `tb_profiler` | nf-core `tbprofiler/profile` | Предсказание лекарственной устойчивости |

#### 7. REPORTS — Отчёты

**Модули:** `multiqc`, `final_table`, `tb_platform_tables`

Этап собирает результаты **всех предыдущих этапов** в единые отчёты:

- **MultiQC** — агрегирует QC-метрики (FastQC, Picard WGS/alignment, samtools stats/flagstat, bcftools stats) в один интерактивный HTML-дашборд
- **FINAL_TABLE.xlsx** — итоговая Excel-таблица со всеми результатами:
  - Образец, линия, споликотип
  - Лекарственная устойчивость
  - Метрики качества
- **TB_PLATFORM_TABLES** — генерирует набор отдельных таблиц для интеграции с платформами:
  - `general.tsv` — общая таблица метрик напрямую из raw Picard / samtools / bcftools outputs, без зависимости от MultiQC
  - `filter.tbmix.tsv` — результаты TB-Mix с аннотацией линий из TBLG
  - `dr.xlsx` — таблица лекарственной устойчивости из TB-Profiler
  - `dr_other_variants.xlsx` — дополнительные варианты, не связанные с устойчивостью
  - `spotyping.total.tsv` — споликотипы всех образцов
  - `rd.tsv` — результаты анализа регионов делеций

#### 8. ANN_TABLE — Когортная таблица аннотаций

**Модули:** nf-core `bcftools/merge`, nf-core `snpeff/snpeff`, локальные `make_table`, `post_process_table`

> **Условие запуска:** ANN_TABLE запускается **только при наличии > 1 образца** и если не задан `--skip_snp_matrix`.

Этап строит **общую матрицу вариантов** для всей когорты:

1. **BCFTOOLS_MERGE** — объединение VCF-файлов всех образцов в один многосемпловый VCF
2. **SNPEFF** — аннотация объединённого VCF через SnpEff
3. **MAKE_TABLE** — извлечение полей через SnpSift extractFields → TSV-таблица с позициями, аллелями, генами и аннотациями
4. **POST_PROCESS_TABLE** — постобработка: добавление информации о генах, стрендах, функциональных категориях; дедупликация записей

---

## Выходные данные

Все результаты сохраняются в `results/`:

| Директория/файл | Содержимое |
|---|---|
| `fastp/` | Отчёты тримминга (HTML, JSON) |
| `fastqc/` | Отчёты качества ридов (HTML, ZIP) |
| `mapped/` | BAM-файлы + метрики дедупликации |
| `stats/picard/` | Picard WGS и alignment метрики |
| `stats/samtools/` | Samtools stats, flagstat |
| `stats/mosdepth/` | Распределение глубины покрытия |
| `vcf/` | VCF-файлы вариантов |
| `annotate_vcf/` | Аннотированные VCF (SnpEff) |
| `snpeff/` | Таблицы аннотаций |
| `spotyping/` | Результаты споликотипирования |
| `lineage/` | Результаты TBLG |
| `rd/` | Результаты RD-анализа |
| `is6110/` | Результаты ISMapper |
| `tb-profiler/` | Результаты TB-Profiler (JSON, TXT) |
| `tb-mix/` | Детекция смешанных инфекций (TSV) |
| `multiqc/` | Агрегированный HTML-отчёт |
| `tb-platform/` | Таблицы для интеграции с платформами (general.tsv, dr.xlsx, rd.tsv и др.) |
| **`FINAL_TABLE.xlsx`** | **Итоговая таблица со всеми результатами** |
| `bad_reads_low_coverage.txt` | Список отфильтрованных образцов |

Для работы SnpEff используется каталог данных, задаваемый параметром `--snpeff_data_dir`.

---

## Инструменты и версии

| Этап | Инструмент | Версия | Контейнер |
|---|---|---|---|
| Тримминг | fastp | 1.1.0 | `staphb/fastp:1.1.0` |
| Тримминг | BBMap repair.sh | 38.96 | bbmap_38.96.sif |
| QC | FastQC | 0.12.1 | `staphb/fastqc:0.12.1` |
| Картирование | BWA + Picard | 3.4.0 | bwa-picard.sif |
| Метрики | Picard | 3.4.0 | `broadinstitute/picard:3.4.0` |
| Варианты | Freebayes | 1.3.10 | freebayes.sif |
| Варианты | BCFtools | 1.22 | `staphb/bcftools:1.22` |
| Варианты | GATK | — | gatk.sif |
| Аннотация | SnpEff | 5.2 | snpeff.sif |
| Генотипирование | SpoTyping | 2.0 | spotyping.sif |
| Генотипирование | ISMapper | — | ismap.sif |
| Генотипирование | Mosdepth | 0.3.11 | mosdpt.sif |
| Генотипирование | TB-Profiler | 6.6.6 | `staphb/tbprofiler:6.6.6` |
| Генотипирование | TBLG | — | tblg.sif |
| RD-анализ | rd.py | — | rd.sif |
| Отчёты | MultiQC | 1.30 | `multiqc/multiqc:v1.30` |

---

## Конфигурация

### Основные параметры (`nextflow.config`)

```groovy
params {
    input              = "run.csv"                    // FASTQ samplesheet
    sra_ids            = null                         // Альтернативный вход: список SRA ID
    outdir             = "./results"                 // Директория результатов
    reference          = "assets/h37rv.fa"           // Референсный геном
    snpeff_data_dir    = "assets/SNPEFF_ANNOTATION/data" // Каталог данных SnpEff
    kraken2_db         = null                        // Первая база Kraken2
    kraken2_db_label   = null                        // Необязательная метка для первой базы
    kraken2_db_2       = null                        // Необязательная вторая база Kraken2
    kraken2_db_label_2 = null                        // Необязательная метка для второй базы
    skip_kraken        = false
    skip_multiqc       = false
    skip_final_reports = false
    mode               = 'link'                      // link | copy

    // Пороги фильтрации
    min_median         = 30
    min_align_pct      = 90
}
```

### Профили запуска

Стандартные runtime-профили:

- `docker` — основной профиль; фактически это режим по умолчанию
- `singularity` — запуск через Singularity / Apptainer
- `conda` — запуск через conda

Для Kubernetes используется отдельный конфиг:

- `-profile docker -c k8s.config`
- `k8s.config` меняет только executor и Kubernetes-specific настройки

### Ресурсы (локальный запуск)

| Лейбл | CPU | Параллельных задач | Память |
|---|---|---|---|
| `solo_cpu` | 1 | 12 | — |
| `multi_cpu` | 4 | 6 | — |
| `small_mem` | — | — | 4 GB |
| `medium_mem` | — | — | 6 GB |
| `big_mem` | — | — | 8 GB |

### Ресурсы (через `-c k8s.config`)

| Лейбл | CPU | Память |
|---|---|---|
| по умолчанию | 3 | задаётся кластером / профилем |

---

## Запуск

### Запуск по FASTQ

```bash
# 1. Подготовить samplesheet
python make_samplesheet.py -i data -o run.csv

# 2. Запустить пайплайн
nextflow run main.nf --input run.csv -resume
```

### Запуск по SRA

```bash
printf "SRR32010433\nERR15166664\n" > sra_ids.txt
nextflow run main.nf --sra_ids sra_ids.txt -resume
```

### Запуск с Kraken2/Bracken

```bash
nextflow run main.nf \
  --input run.csv \
  --kraken2_db /path/to/kraken_db \
  -profile docker -resume
```

### Запуск с двумя базами Kraken2/Bracken

```bash
nextflow run main.nf \
  --input run.csv \
  --kraken2_db /path/to/db1 \
  --kraken2_db_label bacteria \
  --kraken2_db_2 /path/to/db2 \
  --kraken2_db_label_2 viruses \
  -profile docker -resume
```

### Запуск через Singularity

```bash
nextflow run main.nf --input run.csv -profile singularity -resume
```

### Запуск в Kubernetes

```bash
nextflow run main.nf --input run.csv -profile docker -c k8s.config -resume
```

---

## Контейнеры

### Текущая архитектура (Singularity)

Для локальных кастомных модулей доступны отдельные Singularity-контейнеры (.sif):
- **6 кастомных** — собираются из `.def` файлов в `containers/def/`
- **5 публичных** — при необходимости подтягиваются отдельно (fastp, fastqc, bcftools, multiqc, picard)

Скрипт сборки: `build-containers.sh`

Путь к контейнерам задаётся:
- Переменной `NXF_SINGULARITY_CACHEDIR`
- Или по умолчанию: `${baseDir}/containers`

### Кастомные контейнеры

| Контейнер (.sif) | Базовый образ | Инструменты |
|---|---|---|
| build_table.sif | python:3.12-slim | pandas, xlsxwriter, build_final_table.py, build_metrics_table.py |
| tb_platform_tables.sif | python:3.12-slim | profiler_parser.py, deletions_to_csv.py, build_metrics_table.py |
| rd.sif | python:3.11-slim | numpy, pandas, rd.py |
| spotyping.sif | debian:bookworm-slim | SpoTyping 2.0, BLAST 2.12 |
| tblg.sif | python:3.11-slim | tblg (pip) |
| tb_mix.sif | python:3.12-slim | pandas, pysam, tb_mix.py |

---

## Singularity vs Docker: когда что использовать

### Singularity/Apptainer — для HPC-кластеров

Singularity (Apptainer) разработан для HPC-среды:
- Не требует root-прав для запуска контейнеров
- Работает с SLURM, PBS, SGE
- Образы — монолитные .sif файлы (легко переносить)
- Интеграция с общими файловыми системами (NFS, Lustre)

**Используйте Singularity, когда:**
- Запуск на HPC-кластере
- Нет Docker-демона на узлах
- Нужна совместимость с планировщиками задач

### Docker — для серверов и веб-сервисов

Docker — стандарт для серверной контейнеризации:
- Нативная поддержка в Nextflow
- Не нужен `--privileged`
- Образы кэшируются Docker-демоном
- Лёгкий деплой через Docker Compose / Kubernetes

**Используйте Docker, когда:**
- Запуск на сервере с Docker
- Деплой как веб-сервис
- Docker уже установлен и доступен

### Проблема текущей архитектуры

Сейчас пайплайн запускается через Docker, но внутри Docker-контейнера использует Apptainer для запуска .sif файлов:

```
Docker-контейнер (tb-lite:2025-07)
  └── Apptainer запускает .sif файлы
      └── fastp, bwa, freebayes, ...
```

Это **двойная вложенность контейнеров** с проблемами:
1. Требуется `--privileged` (отключает изоляцию Docker)
2. Лишний слой абстракции (оверхед I/O)
3. Docker-образ 8+ ГБ (все .sif внутри)
4. Архитектурная избыточность

### Рекомендация

Для сервера: перейти на нативный Docker-режим Nextflow. Это устранит все перечисленные проблемы. Singularity оставить для HPC.
