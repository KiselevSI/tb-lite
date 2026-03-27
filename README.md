# TB-Lite: Nextflow-пайплайн для геномного анализа M. tuberculosis

## Обзор

TB-Lite — это Nextflow-пайплайн для полного геномного анализа *Mycobacterium tuberculosis* (палочка Коха). Пайплайн принимает FASTQ-файлы секвенирования и выполняет полный цикл: от тримминга ридов до определения линии, споликотипа и лекарственной устойчивости.

**Референсный геном:** H37Rv (*M. tuberculosis*)

---

## Структура проекта

```
tb-lite/
├── main.nf                    # Главный workflow (точка входа)
├── nextflow.config            # Конфигурация: параметры, профили, контейнеры
├── Dockerfile                 # Docker-образ для деплоя на сервере
├── build-containers.sh        # Скрипт сборки Singularity-контейнеров
├── run.csv                    # Входной файл с образцами (пример)
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
├── modules/                   # 26 процессов (атомарные единицы)
│   ├── trimming/              # fastp, repair_reads
│   ├── qc/                    # fastqc
│   ├── mapping/               # bwa_index, bwa_picard
│   ├── call_variants/         # freebayes, bcftools, gatk, snpeff, rename_chr, bcftools_stats
│   ├── filter/                # tb_mix, map_stats, samtools_stats
│   ├── genotyping/            # spotyping, ismapper, mosdepth, rd, tblg, tb_profiler
│   ├── reports/               # multiqc, final_table
│   └── ann_table/             # merge, snpeff_ann, make_table, post_process_table
│
├── assets/                    # Референсные данные
│   ├── h37rv.fa               # Референсный геном H37Rv
│   ├── rd/RD.bed              # Регионы различий (RD)
│   ├── ismap/is6110.fasta     # Последовательность IS6110
│   ├── ismap/h37rv.gbk        # GenBank-файл H37Rv
│   ├── chr_name/chr.txt       # Маппинг имён хромосом
│   ├── tbmix/levels.tsv       # Уровни для детекции смешанных инфекций
│   ├── vcf2table/feature_table.tsv
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

### Файл образцов (`run.csv`)

CSV-файл с 4 колонками:

| Колонка  | Описание |
|----------|----------|
| `Sample` | Имя образца (уникальное) |
| `R1`     | Полный путь к FASTQ R1 |
| `R2`     | Полный путь к FASTQ R2 (пусто для single-end) |
| `Layout` | `PAIRED` или `SOLO` |

**Пример:**
```csv
Sample,R1,R2,Layout
ERR123,/data/ERR123_1.fastq.gz,/data/ERR123_2.fastq.gz,PAIRED
SRR456,/data/SRR456.fastq.gz,,SOLO
```

---

## Конвейер обработки

### Общая схема (8 этапов)

```
run.csv
    │
    ▼
┌─────────────────────────────────────┐
│ 1. TRIMMING (fastp)                 │  Тримминг адаптеров, polyG, качества
│    Модули: fastp, repair_reads      │  repair_reads — опционально (BBMap)
└─────────────┬───────────────────────┘
              │
    ┌─────────┴─────────┐
    ▼                   ▼
┌──────────┐    ┌──────────────────────┐
│ 2. QC    │    │ 3. MAPPING           │
│ (FastQC) │    │ BWA-mem → SAMtools   │
│          │    │ → Picard MarkDups    │
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

### Подробное описание этапов

#### 1. TRIMMING — Тримминг ридов

**Модули:** `fastp`, `repair_reads` (опционально)

- **fastp** — обрезка адаптеров (автодетекция), polyG-хвостов, фильтрация по качеству
- **repair_reads** (BBMap repair.sh) — восстановление сломанных пар (если нужно)

**Выход:** Триммированные FASTQ-файлы

#### 2. QC — Контроль качества

**Модуль:** `fastqc`

- Генерация отчётов FastQC для каждого образца
- Оценка качества секвенирования, GC-состава, длин ридов

#### 3. MAPPING — Картирование

**Модули:** `bwa_index`, `bwa_picard`

- **BWA-mem** — выравнивание ридов на референс H37Rv
- **SAMtools** — сортировка и индексация BAM
- **Picard MarkDuplicates** — маркировка PCR-дупликатов

**Выход:** Отсортированные, дедуплицированные BAM-файлы

#### 4. FILTER — Фильтрация образцов

**Модули:** `tb_mix`, `map_stats`, `samtools_stats`

- **TB-Mix** (tb_mix.py) — детекция смешанных инфекций (несколько штаммов в образце)
- **Picard CollectWgsMetrics** — метрики покрытия генома
- **Picard CollectAlignmentSummaryMetrics** — метрики выравнивания
- **Samtools stats/flagstat** — статистика BAM-файлов

**Фильтры (настраиваемые):**
- Медианное покрытие >= 10X (`min_median`) — **активный фильтр**
- Среднее покрытие >= 10X (`min_mean_cov`)
- % выравненных >= 80% (`min_align_pct`)

Образцы, не прошедшие фильтр → `bad_reads_low_coverage.txt`

#### 5. VARIANT CALLING — Вызов вариантов

**Модули:** `freebayes` (основной), `bcftools`, `gatk` (альтернативные), `snpeff`, `rename_chr`, `bcftools_stats`

**Freebayes** (основной вариант-каллер):
- Гаплоидный режим (`--ploidy 1`) — TB не имеет диплоидного генома
- Минимальное качество картирования: 20
- Минимальное качество баз: 20
- Минимальное покрытие: 10X
- Минимальная частота альтернативного аллеля: 80%

**Постобработка:**
- BCFtools norm — нормализация VCF
- BCFtools filter — фильтрация по QUAL >= 10
- SnpEff — аннотация вариантов (эффект на гены, аминокислотные замены)
- Переименование хромосом (chr.txt) для совместимости

#### 6. GENOTYPING — Генотипирование

**Модули** (запускаются параллельно):

| Модуль | Инструмент | Что делает |
|--------|-----------|------------|
| `spotyping` | SpoTyping | Определяет споликотип (24-спейсерный бинарный код) |
| `ismapper` | ISMapper | Картирование IS6110 (только для парных ридов) |
| `mosdepth` | Mosdepth | Быстрый расчёт глубины покрытия |
| `rd` | rd.py | Анализ Region of Difference (RD) для линии |
| `tblg` | TBLG | Классификация линии (L1–L7) по VCF |
| `tb_profiler` | TB-Profiler | Предсказание лекарственной устойчивости |

#### 7. REPORTS — Отчёты

**Модули:** `multiqc`, `final_table`

- **MultiQC** — агрегирует все отчёты QC в один HTML-дашборд
- **FINAL_TABLE.xlsx** — итоговая Excel-таблица со всеми результатами:
  - Образец, линия, споликотип
  - Лекарственная устойчивость
  - Метрики качества

#### 8. ANN_TABLE — Когортная таблица аннотаций

**Модули:** `merge`, `snpeff_ann`, `make_table`, `post_process_table`

- Объединение VCF всех образцов когорты
- Аннотация объединённого VCF через SnpEff
- Извлечение полей через SnpSift → TSV-таблица
- Постобработка: добавление генов, стрендов, дедупликация

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
| **`FINAL_TABLE.xlsx`** | **Итоговая таблица со всеми результатами** |
| `bad_reads_low_coverage.txt` | Список отфильтрованных образцов |

---

## Инструменты и версии

| Этап | Инструмент | Версия | Контейнер |
|---|---|---|---|
| Тримминг | fastp | 0.24.1 | `staphb/fastp:0.24.1` |
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
| Генотипирование | TB-Profiler | 6.6.5 | tbp.sif |
| Генотипирование | TBLG | — | tblg.sif |
| RD-анализ | rd.py | — | rd.sif |
| Отчёты | MultiQC | 1.30 | `multiqc/multiqc:v1.30` |

---

## Конфигурация

### Основные параметры (`nextflow.config`)

```groovy
params {
    outdir        = "./results"           // Директория результатов
    samples       = "run.csv"             // Файл с образцами
    reference     = "assets/h37rv.fa"     // Референсный геном
    mode          = 'link'                // link | copy (режим создания файлов)

    // Пороги фильтрации
    min_median    = 10    // Мин. медианное покрытие (активный)
    min_mean_cov  = 10    // Мин. среднее покрытие
    min_align_pct = 80    // Мин. % выравненных ридов
}
```

### Ресурсы (профиль `local`)

| Лейбл | CPU | Параллельных задач | Память |
|---|---|---|---|
| `solo_cpu` | 1 | 12 | — |
| `multi_cpu` | 2 | 6 | — |
| `small_mem` | — | — | 4 GB |
| `medium_mem` | — | — | 6 GB |
| `big_mem` | — | — | 8 GB |

---

## Запуск

### Локальный запуск (Singularity)

```bash
# 1. Собрать контейнеры (один раз)
bash build-containers.sh

# 2. Подготовить run.csv с путями к FASTQ-файлам

# 3. Запустить пайплайн
nextflow run main.nf -profile local -resume
```

### Запуск через Docker (на сервере)

```bash
docker run --rm -it --privileged \
    -v $(pwd)/data:/workspace/data \
    -v $(pwd)/results:/workspace/results \
    tb-lite:2025-07 \
    run main.nf -with-singularity -resume
```

> **Важно:** Флаг `--privileged` необходим для запуска Apptainer внутри Docker.
> Рекомендуется перейти на нативный Docker-режим (см. раздел ниже).

---

## Контейнеры

### Текущая архитектура (Singularity)

Каждый процесс запускается в своём Singularity-контейнере (.sif):
- **13 кастомных** — собираются из `.def` файлов в `containers/def/`
- **5 публичных** — скачиваются из Docker Hub (fastp, fastqc, bcftools, multiqc, picard)

Скрипт сборки: `build-containers.sh`

Путь к контейнерам задаётся:
- Переменной `NXF_SINGULARITY_CACHEDIR`
- Или по умолчанию: `${baseDir}/containers`

### Кастомные контейнеры

| Контейнер (.sif) | Базовый образ | Инструменты |
|---|---|---|
| ann_table.sif | python:3.12-slim | pandas, xlsxwriter, кастомные скрипты |
| build_table.sif | python:3.12-slim | pandas, xlsxwriter, build_final_table.py |
| bwa-picard.sif | ubuntu:20.04 | bwa, samtools, Picard 3.4.0 |
| freebayes.sif | ubuntu:24.04 | freebayes 1.3.10, bcftools |
| ismap.sif | python:3.8-slim | bwa, samtools, bedtools, blast, IS_mapper |
| mosdpt.sif | ubuntu:22.04 | mosdepth 0.3.11 |
| rd.sif | python:3.11-slim | numpy, pandas, rd.py |
| snpeff.sif | miniconda3 | snpEff 5.2 + БД H37Rv |
| snpeff_table.sif | ubuntu:24.04 | snpEff, SnpSift, bcftools |
| spotyping.sif | debian:bookworm-slim | SpoTyping 2.0, BLAST 2.12 |
| tblg.sif | python:3.11-slim | tblg (pip) |
| tb-lineage.sif | python:3.12-slim | pandas, pysam, tb-lineage.py |
| tb_mix.sif | python:3.12-slim | pandas, pysam, tb_mix.py |
| tbp.sif | miniconda3 | TB-Profiler 6.6.5 |
| vcf2table.sif | python:3.12-slim | pandas, pysam, vcf2table.py |

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
