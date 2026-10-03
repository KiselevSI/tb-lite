# TB-Lite

Nextflow DSL2-пайплайн для WGS-анализа *Mycobacterium tuberculosis*: от FASTQ
(или SRA accession) до лекарственной устойчивости, линий, сполиготипа, RD,
IS6110, SNP-матрицы и таблиц для TB Platform. Референс — H37Rv.

Этот файл — быстрый старт: как запускать и что получится на выходе.
Как пайплайн устроен внутри, какие инструменты и параметры используются —
в [DOCUMENTATION.md](DOCUMENTATION.md).

## Содержание

- [Требования](#требования)
- [Быстрый старт](#быстрый-старт)
- [Входные данные](#входные-данные)
- [Режимы запуска](#режимы-запуска)
  - [1. FASTQ](#1-fastq)
  - [2. SRA](#2-sra)
  - [3. Kraken2/Bracken](#3-kraken2bracken)
  - [4. Batch-режим (большие наборы)](#4-batch-режим-большие-наборы)
  - [5. Сборка отчётов по готовому каталогу](#5-сборка-отчётов-по-готовому-каталогу)
  - [6. Только аннотация VCF](#6-только-аннотация-vcf)
  - [7. SNP-матрица из готовых VCF](#7-snp-матрица-из-готовых-vcf)
- [Профили запуска](#профили-запуска)
- [Что получается на выходе](#что-получается-на-выходе)
- [Основные параметры](#основные-параметры)
- [Перезапуск и частые проблемы](#перезапуск-и-частые-проблемы)

## Требования

- Nextflow ≥ 23.04 (проверено на 25.10) и Java 17+.
- Один из runtime: Docker, Singularity/Apptainer или Conda.
- Для Docker нужно один раз собрать локальные образы `tb-lite/*` (в публичных
  реестрах их нет), а после обновления пайплайна пересобрать их:

  ```bash
  bash build-docker-images.sh
  ```

  Если образ не собран, Nextflow падает с
  `pull access denied for tb-lite/...` — см. [частые проблемы](#перезапуск-и-частые-проблемы).
  Про Singularity и локальные образы см.
  [DOCUMENTATION.md](DOCUMENTATION.md#13-известные-ограничения).

## Быстрый старт

```bash
# samplesheet из папки с *.fastq.gz / *.fq.gz
python make_samplesheet.py -i /data/fastq -o run.csv

nextflow run main.nf -profile docker \
  --input run.csv \
  --outdir results \
  -resume
```

Главный результат — `results/Reports/general/FINAL_TABLE.xlsx`.

## Входные данные

| Режим | Параметр | Формат |
|---|---|---|
| FASTQ | `--input` | CSV `sample,fastq_1,fastq_2` |
| SRA | `--sra_ids` | TXT, один accession на строку |
| Аннотация VCF | `--vcf_list` | CSV `sample,vcf` |
| SNP-матрица из VCF | `--vcf_input` | CSV `sample,vcf` |

**FASTQ samplesheet.** Нужны полные пути, только gzip (`*.fastq.gz` или
`*.fq.gz`). Для single-end колонка `fastq_2` остаётся пустой. Старый формат
`Sample,R1,R2,Layout` тоже принимается.

```csv
sample,fastq_1,fastq_2
ERR123,/data/ERR123_1.fastq.gz,/data/ERR123_2.fastq.gz
SRR456,/data/SRR456.fastq.gz,
```

**SRA.** Пустые строки и строки с `#` пропускаются, дубликаты убираются.

```text
SRR32010433
ERR15166664
```

**VCF-список.** В каждом VCF должна быть ровно одна sample-колонка. Имя
образца — только буквы, цифры, `.`, `_`, `-`.

```csv
sample,vcf
sample1,/data/sample1.vcf
sample2,/data/sample2.vcf.gz
```

## Режимы запуска

### 1. FASTQ

```bash
nextflow run main.nf -profile docker \
  --input run.csv \
  --outdir results \
  -resume
```

Полный анализ: тримминг → QC → картирование → фильтр по покрытию →
варианты → генотипирование → отчёты. Если в наборе больше одного образца,
строится ещё и когортная SNP-матрица.

### 2. SRA

```bash
nextflow run main.nf -profile docker \
  --sra_ids sra_ids.txt \
  --outdir results \
  -resume
```

То же, что режим FASTQ, но риды сначала скачиваются из SRA
(`prefetch` + `fasterq-dump`). Accession, у которых после скачивания не 1 и не 2
FASTQ-файла, пропускаются и попадают в `Reports/general/unsupported_sra_layout.txt`.

### 3. Kraken2/Bracken

По умолчанию Kraken не запускается: он включается, если указать `--kraken2_db`.
Баз может быть одна или две.

```bash
nextflow run main.nf -profile docker \
  --input run.csv \
  --kraken2_db   /db/kraken_standard --kraken2_db_label   ALL \
  --kraken2_db_2 /db/kraken_myco     --kraken2_db_label_2 ONLY_MYCOBACTERIUM \
  --outdir results \
  -resume
```

Если метка не задана, она берётся из имени каталога базы. Метки двух баз
не должны совпадать.

### 4. Batch-режим (большие наборы)

Для сотен и тысяч образцов. `run_batches.sh` делит вход на батчи, прогоняет
их по очереди в один общий `outdir`, после каждого успешного батча чистит
`work/` и в конце собирает общие отчёты.

```bash
./run_batches.sh \
  --pipeline /home/zerg/git/tb-lite \
  --input /data/run.csv \
  --batch-size 500 \
  --profile docker \
  --outdir /data/results_batches \
  --workdir /data/work_batches
```

| Опция | По умолчанию | Назначение |
|---|---|---|
| `--pipeline` | — (обязательна) | Каталог с `main.nf` |
| `--input` | — (обязательна) | Samplesheet CSV или список SRA |
| `--input-mode` | `auto` | `fastq`/`sra`; `auto`: если в первой строке есть запятая — FASTQ |
| `--batch-size` | `500` | Образцов в батче |
| `--profile` | **`conda`** | `docker`, `singularity`, `conda`, `local` |
| `--outdir` | `./results` | Общий каталог результатов |
| `--workdir` | `./work` | Work-каталог Nextflow (**очищается после каждого батча**) |
| `--resume-from N` | `1` | Продолжить с батча N |
| `--with-kraken`, `--kraken2_db*` | выкл. | Kraken, как в режиме 3 |
| `--benchmark` | выкл. | Трассировка Nextflow и сводка производительности в `<outdir>/benchmark/` |

Разбивка сохраняется в `./.batches/`, журнал батчей — в `./batches.log`
(оба в текущем каталоге). Если батч упал, скрипт печатает готовую команду
с `--resume-from` для продолжения, а `work/` не трогает.

> Не используйте `--mode link`/`symlink` в batch-режиме: `work/` удаляется
> после каждого батча, и ссылки в `outdir` станут битыми.

С `--benchmark` для каждого батча сохраняются `trace.tsv`, `report.html` и
`timeline.html` Nextflow, а в конце считается сводка: образцов в час/сутки,
CPU-часы и пиковая память на образец, время по процессам, конфигурация
машины. Главный файл — `<outdir>/benchmark/summary.md`. Подробности — в
[DOCUMENTATION.md](DOCUMENTATION.md#бенчмарк---benchmark).

### 5. Сборка отчётов по готовому каталогу

`batch_reports.nf` не запускает анализ, а собирает отчёты из того, что уже
лежит в `outdir`. `run_batches.sh` вызывает его сам. Вручную он нужен, чтобы
пересобрать отчёты или получить SNP-таблицы для TB Platform после обычного
запуска.

```bash
# Все отчёты + SNP-таблицы
nextflow run batch_reports.nf -profile docker --outdir /data/results_batches

# Только SNP-таблицы для TB Platform
nextflow run batch_reports.nf -entry SNP_DB_ONLY -profile docker --outdir /data/results
```

Полная сборка падает, если в `outdir` нет хотя бы одного из ожидаемых типов
файлов (например, `fastqc/` при запуске с `--skip_qc`). В таком случае
используйте `-entry SNP_DB_ONLY`. Для Kraken передайте те же `--kraken2_db*`,
иначе добавьте `--skip_kraken`.

### 6. Только аннотация VCF

SnpEff-аннотация готовых per-sample VCF без WGS-анализа:

```bash
nextflow run main.nf -profile docker \
  --vcf_annotation_only \
  --vcf_list vcf_samples.csv \
  --outdir results_vcf_annotation
```

`--input`/`--sra_ids` в этом режиме не задаются.

### 7. SNP-матрица из готовых VCF

Когортная SNP-матрица из готовых VCF (аннотировать их заранее не нужно):

```bash
nextflow run snp_matrix.nf -profile docker \
  --vcf_input vcf_samples.csv \
  --outdir results_snp_matrix
```

CSV можно собрать из каталогов с VCF:

```bash
# sample ID из имени файла — быстро
python make_snp_matrix_csv.py -i /data/vcf1 /data/vcf2 -o vcf_samples.csv \
  --sample-source filename --jobs 16

# sample ID из заголовка VCF — надёжнее, но открывает каждый файл
python make_snp_matrix_csv.py -i /data/vcf -o vcf_samples.csv \
  --sample-source header --jobs 16
```

В матрицу не попадают образцы с 5000 и более вариантами (обычно это
контаминация или не MTBC). Если после этого осталось меньше двух образцов,
матрица не строится.

## Профили запуска

| Профиль | Когда |
|---|---|
| `-profile docker` | Обычный сервер с Docker |
| `-profile singularity` | HPC/кластер с Apptainer/Singularity |
| `-profile conda` | Нет контейнеров; окружения из `environment.yml` модулей |
| `-profile docker -c k8s.config` | Kubernetes (executor `k8s`) |

Без `-profile` используется Docker, но лучше указывать профиль явно.

## Что получается на выходе

### Режимы 1–3 (FASTQ, SRA, Kraken)

| Путь | Что это |
|---|---|
| `Reports/general/FINAL_TABLE.xlsx` | **Главная сводная таблица**: все образцы после fastp, QC-метрики, TB-Mix, линия, сполиготип, устойчивость, Kraken top-5 |
| `Reports/general/drug_resist.xlsx` | Лекарственная устойчивость (TB-Profiler) |
| `Reports/general/bad_reads_low_coverage.txt` | Отбракованные по % картирования / медиане покрытия |
| `Reports/general/bad_reads_invalid_fastq.txt` | Отбракованные из-за пустых или битых FASTQ |
| `Reports/general/unsupported_sra_layout.txt` | Только SRA: accession с неподдерживаемым layout |
| `Reports/tb-platform/` | Таблицы для загрузки в TB Platform (только прошедшие фильтр образцы) |
| `Reports/tb-platform/is6110/` | Нормализованные таблицы вставок IS6110 |
| `Reports/snp_matrix/FINAL_ANNOTATION_TABLE.tsv` | Когортная SNP-матрица (если образцов ≥ 2) |
| `multiqc/TB-Lite-QC_multiqc_report.html` | Сводный QC-отчёт MultiQC |
| `versions.txt` | Версии всех использованных программ |
| `snp_db/shards/` | Заготовки SNP-таблиц TB Platform; готовые таблицы — через режим 5 |

Кроме того, в `outdir` лежат результаты по каждому образцу: `fastp/`,
`fastqc/`, `stats/`, `vcf/`, `annotate_vcf/`, `tb-mix/`, `lineage/`,
`spotyping/`, `rd/`, `tb-profiler/`, `is6110/`, `kraken2/` (если включён Kraken).
Полный список — в [DOCUMENTATION.md](DOCUMENTATION.md#11-выходные-каталоги).
BAM-файлы не сохраняются.

### Режим 4 (batch)

Каталоги по образцам — те же, что в режимах 1–3, но общие для всех батчей.
После последнего батча собираются:

| Путь | Что это |
|---|---|
| `Reports/general/`, `Reports/tb-platform/` | Общие отчёты по всем батчам |
| `Reports/tb-platform/snp/` | Готовые SNP-таблицы и `import_snp.sql` для TB Platform |
| `Reports/general/bad_reads_*.txt` | Объединённые списки отбракованных |
| `batch_reports/filter/` | Те же списки отбракованных по каждому батчу |
| `benchmark/` | Только с `--benchmark`: trace/report/timeline по батчам и сводка производительности |

В batch-режиме **не строятся** MultiQC и когортная SNP-матрица. Матрицу
при необходимости строят отдельно (режим 7) по `vcf/`.

### Режим 5 (`batch_reports.nf`)

Всё из `Reports/` как в batch-режиме. С `-entry SNP_DB_ONLY` — только
`Reports/tb-platform/snp/`.

### Режим 6 (аннотация VCF)

| Путь | Что это |
|---|---|
| `annotate_vcf/<sample>/<sample>.annotated.ann.vcf` | VCF с SnpEff-аннотацией, sample-колонка переименована по CSV |
| `annotate_vcf/<sample>/snpEff_summary.html` и др. | Отчёты SnpEff |
| `versions.txt` | Версии программ |

### Режим 7 (SNP-матрица)

| Путь | Что это |
|---|---|
| `Reports/snp_matrix/FINAL_ANNOTATION_TABLE.tsv` | SNP-матрица с аннотацией |
| `Reports/snp_matrix/cohort.*` | Промежуточные когортные VCF |

Как загрузить результаты в базу TB Platform, описано в
[DOCUMENTATION.md](DOCUMENTATION.md#7-загрузка-в-tb-platform).

## Основные параметры

| Параметр | По умолчанию | Назначение |
|---|---|---|
| `--outdir` | `./results2` | Каталог результатов |
| `--min_align_pct` | `90` | Минимальный % картированных ридов |
| `--min_median` | `30` | Минимальная медиана покрытия |
| `--mode` | `copy` | Как публиковать файлы (`copy`, `link`, `symlink`) |
| `--skip_qc` | `false` | Без FastQC |
| `--skip_kraken` | `false` | Без Kraken, даже если задана база |
| `--skip_multiqc` | `false` | Без MultiQC |
| `--skip_final_reports` | `false` | Без `Reports/general` и `Reports/tb-platform` |
| `--skip_snp_matrix` | `false` | Без когортной SNP-матрицы |
| `--skip_snp_db` | `false` | Без SNP-таблиц TB Platform |

Полный список (референс, SnpEff, базы, ресурсы) — в
[DOCUMENTATION.md](DOCUMENTATION.md#9-параметры). Справка: `nextflow run main.nf --help`.

## Перезапуск и частые проблемы

- **Перезапуск после падения** — та же команда с `-resume`: готовые задачи
  возьмутся из кэша. В batch-режиме — `--resume-from N`.
- **`pull access denied for tb-lite/<имя>`** — на машине нет локального
  образа нужной версии. Соберите его: `bash build-docker-images.sh`.
- **Образца нет в отчётах** — проверьте `Reports/general/bad_reads_*.txt`.
  Отбракованные образцы остаются в `FINAL_TABLE.xlsx` с пустыми полями
  после фильтра, но в `Reports/tb-platform/` не попадают.
- **ISMapper/IS6110 пустые** — ISMapper работает только с paired-end ридами.
