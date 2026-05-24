# TB-Lite: Nextflow-пайплайн для геномного анализа *M. tuberculosis*

## Обзор

TB-Lite — это Nextflow DSL2-пайплайн для полного WGS-анализа *Mycobacterium tuberculosis*. Пайплайн принимает либо локальный samplesheet с `FASTQ.gz`, либо список SRA accession ID и выполняет полный цикл анализа: QC, картирование, фильтрацию образцов, вызов вариантов, генотипирование, предсказание лекарственной устойчивости, опциональную таксономическую классификацию Kraken2/Bracken и сборку итоговых отчётов.

Референс по умолчанию: H37Rv.

## Ключевые особенности

- DSL2-пайплайн на базе nf-core modules и локальных TB-специфичных модулей.
- Поддержка двух режимов входа: `--input` для локальных FASTQ и `--sra_ids` для SRA.
- Стандартные runtime-профили: `docker`, `singularity`, `conda`.
- Опциональный Kraken2/Bracken с одной или двумя базами.
- Автоматическая фильтрация образцов по `% mapped` и `median coverage`.
- Итоговый `FINAL_TABLE.xlsx` включает все образцы, дошедшие до `fastp`; для отфильтрованных downstream-поля остаются пустыми.
- Когортная SNP-матрица для наборов из более чем одного образца.
- Batch-режим через `run_batches.sh` с финальной агрегацией одного общего `Reports/` и одного общего `multiqc/`.

## Структура проекта

```text
tb-lite/
├── main.nf
├── batch_reports.nf
├── nextflow.config
├── conf/
│   ├── base.config
│   └── modules.config
├── workflows/
│   └── tblite.nf
├── subworkflows/
│   └── local/
├── modules/
│   ├── nf-core/
│   └── local/
├── lib/
│   └── WorkflowMain.groovy
├── bin/
├── assets/
│   ├── h37rv.fa
│   ├── h37rv.gbk
│   ├── SNPEFF_ANNOTATION/
│   │   ├── data/
│   │   └── h37rv_feature_table.txt
│   ├── multiqc/
│   ├── rd/
│   ├── ismap/
│   ├── chr_name/
│   └── tbmix/
├── containers/
│   ├── dockerfiles/
│   └── def/
├── run_batches.sh
├── build-docker-images.sh
├── build-containers.sh
└── make_samplesheet.py
```

`containers/dockerfiles/` и `containers/def/` содержат runtime-описания для локальных модулей. Это не означает, что обычный запуск использует "Docker внутри Apptainer": реальный runtime выбирается профилем Nextflow.

## Входные данные

### 1. FASTQ samplesheet

Параметр: `--input`

CSV-файл с колонками:

| Колонка | Описание |
|---|---|
| `sample` | Идентификатор образца |
| `fastq_1` | Полный путь к `*.fastq.gz` или `*.fq.gz` для R1 или single-end |
| `fastq_2` | Полный путь к R2, пусто для single-end |

Совместимость со старым форматом `Sample,R1,R2,Layout` сохранена.
Поддерживаются только gzipped FASTQ: `*.fastq.gz` или `*.fq.gz`.

Пример:

```csv
sample,fastq_1,fastq_2
ERR123,/data/ERR123_1.fastq.gz,/data/ERR123_2.fastq.gz
SRR456,/data/SRR456.fastq.gz,
```

### 2. Список SRA accession

Параметр: `--sra_ids`

Текстовый файл с одним accession ID на строку.

Пример:

```text
SRR32010433
ERR15166664
```

## Логика пайплайна

### Основной поток

1. `TRIMMING`
   `fastp` обрезает адаптеры, polyG и низкокачественные хвосты.
2. `QC`
   `FastQC` строит отчёты по trimmed reads, если не задан `--skip_qc`.
3. `MAPPING`
   `bwa mem` и `Picard MarkDuplicates` строят дедуплицированный BAM.
4. `FILTER`
   `TB-Mix`, `Picard CollectWgsMetrics`, `Picard CollectAlignmentSummaryMetrics`, `samtools stats` и `samtools flagstat` оценивают качество образца.
5. `KRAKEN`
   Опциональная ветка Kraken2/Bracken запускается на ридах сразу после `fastp`, до sample filtering.
6. `CALLVAR`
   Для прошедших фильтр образцов запускаются `Freebayes`, `BCFtools` и `SnpEff`.
7. `GENOTYPE`
   Для прошедших фильтр образцов запускаются `SpoTyping`, `ISMapper`, `Mosdepth`, `RD`, `TBLG`, `TB-Profiler`.
8. `REPORTS`
   Собираются `MultiQC`, `Reports/general`, `Reports/tb-platform` и, при необходимости, `Reports/snp_matrix`.

### Фильтрация образцов

По умолчанию образец считается "хорошим", если выполняются оба условия:

- `reads_mapped_percent >= 90` (`--min_align_pct`)
- `median_coverage >= 30` (`--min_median`)

Непрошедшие образцы попадают в `Reports/general/bad_reads_low_coverage.txt` и не идут в downstream-ветки `CALLVAR`, `GENOTYPE` и `ANN_TABLE`.

При этом:

- Kraken, если включён, всё равно считается для них, потому что запускается раньше фильтрации.
- В `FINAL_TABLE.xlsx` строка для такого образца сохраняется, но downstream-поля будут пустыми.

## Отчёты и как они формируются

### `MultiQC`

`MultiQC` собирается из опубликованных sample-level артефактов:

- `fastp` JSON
- `FastQC` ZIP
- `Picard CollectWgsMetrics`
- `Picard CollectAlignmentSummaryMetrics`
- `samtools stats`
- `samtools flagstat`
- `bcftools stats`
- Kraken2 reports, если Kraken включён

Важно: `FINAL_TABLE.xlsx` не строится из `MultiQC`. `MultiQC` и финальные табличные отчёты — это параллельные отчётные ветки.

### `FINAL_TABLE.xlsx`

`FINAL_TABLE.xlsx` публикуется в `Reports/general/FINAL_TABLE.xlsx` и собирается напрямую из опубликованных raw/published outputs:

- полный список образцов после `fastp`
- `general.tsv` из Picard / samtools / bcftools метрик
- `tbmix.total.tsv`
- `filter.tbmix.tsv`
- `spotyping.total.tsv`
- `tblg.total.tsv`
- `drug_resist.xlsx`
- `kraken.top_hits.tsv`, если Kraken включён

Особенности:

- Базой служит список всех образцов после `fastp`.
- Для образцов, отфильтрованных позже, строки сохраняются.
- При включённом Kraken добавляются колонки `*_top1..top5` для каждой Kraken DB.
- В каждой Kraken-ячейке указывается организм и доля из `*_frac`, округлённая до двух знаков.

### `Reports/tb-platform`

Процесс `TB_PLATFORM_TABLES` публикует отдельные файлы в `Reports/tb-platform/`:

- `general.tsv`
- `filter.tbmix.tsv`
- `drug_resist.xlsx`
- `drug_resist_and_uncertain.xlsx`
- `spotyping.total.tsv`
- `rd.tsv`

### `Reports/snp_matrix`

Когортная матрица строится, если в наборе больше одного образца и не задан `--skip_snp_matrix`.

Основной финальный файл:

- `Reports/snp_matrix/FINAL_ANNOTATION_TABLE.tsv`

Также публикуются промежуточные cohort-level VCF/annotation файлы.

### VCF annotation-only

Если уже есть per-sample VCF и нужно только получить SnpEff-annotated VCF без полного WGS-запуска:

```bash
nextflow run . \
  -profile docker \
  --vcf_annotation_only \
  --vcf_list vcf_samples.csv \
  --outdir results_vcf_annotation
```

Формат `vcf_samples.csv`:

```csv
sample,vcf
sample1,/data/sample1.vcf
sample2,/data/sample2.vcf.gz
```

Каждый входной VCF должен содержать ровно один sample column. Пайплайн переименует sample column в значение из колонки `sample` и опубликует результат в `annotate_vcf/<sample>/`.

### Standalone SNP matrix из VCF

Если уже есть per-sample VCF, SNP-матрицу можно построить отдельным entrypoint без полного WGS-запуска:

```bash
nextflow run snp_matrix.nf \
  -profile docker \
  --vcf_input vcf_samples.csv \
  --outdir results_snp_matrix
```

Минимальный вход — обычный однообразцовый VCF (`.vcf` или `.vcf.gz`) для каждого образца. Per-sample annotated VCF не требуется: workflow сначала объединяет VCF в cohort VCF, затем запускает SnpEff и строит `FINAL_ANNOTATION_TABLE.tsv`.

Формат `vcf_samples.csv`:

```csv
sample,vcf
sample1,/data/sample1.vcf
sample2,/data/sample2.vcf.gz
```

CSV можно создать из директории с VCF:

```bash
python make_snp_matrix_csv.py -i /data/vcf -o vcf_samples.csv
```

Можно передать несколько директорий:

```bash
python make_snp_matrix_csv.py \
  -i ../VCF/ /data6/bio/MolGenMicro/TBGenoPipe/results/VCF2/VCF \
  -o vcf_samples.csv
```

## Выходные данные

### Итоговые отчёты

| Путь | Назначение |
|---|---|
| `Reports/general/FINAL_TABLE.xlsx` | Главный итоговый Excel-отчёт |
| `Reports/general/drug_resist.xlsx` | Drug-resistance таблица из `profiler_parser.py` |
| `Reports/general/bad_reads_low_coverage.txt` | Непрошедшие фильтр образцы |
| `Reports/tb-platform/` | Отдельные платформенные TSV/XLSX-таблицы |
| `Reports/snp_matrix/FINAL_ANNOTATION_TABLE.tsv` | Когортная SNP-матрица |
| `multiqc/multiqc_report.html` | Итоговый HTML-отчёт MultiQC |

### Sample-level published outputs

Эти директории важны не только для отладки, но и для batch aggregation:

| Путь | Назначение |
|---|---|
| `fastp/<sample>/` | `fastp` JSON/HTML |
| `fastqc/<sample>/` | FastQC ZIP/HTML |
| `mapped/<sample>/` | BAM и BAM index после дедупликации |
| `stats/picard/wgs/<sample>/` | `CollectWgsMetrics.coverage_metrics` |
| `stats/picard/alignment/<sample>/` | `CollectAlignmentSummaryMetrics` |
| `stats/samtools/stats/<sample>/` | `samtools stats` |
| `stats/samtools/flagstat/<sample>/` | `samtools flagstat` |
| `stats/bcftools/` | `*.bcftools_stats.txt` |
| `stats/mosdepth/<sample>/` | Mosdepth outputs |
| `vcf/<sample>/` | filtered VCF |
| `annotate_vcf/<sample>/` | SnpEff-annotated VCF |
| `tb-mix/` | результаты TB-Mix |
| `spotyping/<sample>/` | результаты SpoTyping |
| `lineage/` | lineage-таблицы TBLG |
| `rd/<sample>/` | RD tables |
| `tb-profiler/drug-resist/<sample>/results/` | `*.results.json` и другие outputs TB-Profiler |
| `is6110/paired/<sample>/` | ISMapper outputs |
| `kraken2/kraken2/<db_label>/<sample>/` | Kraken2 per-sample outputs |
| `kraken2/bracken/<db_label>/<sample>/` | Bracken per-sample outputs |
| `kraken2/combined/` | combined all-sample tables по каждой Kraken DB |

## Runtime-профили

TB-Lite поддерживает три стандартных профиля:

- `docker`
- `singularity`
- `conda`

Если профиль не указан, в текущем `nextflow.config` Docker runtime остаётся включён по умолчанию. Для воспроизводимых запусков лучше указывать профиль явно.

### Рекомендации

- `-profile docker`
  Обычный серверный запуск с Docker.
- `-profile singularity`
  HPC и кластеры с Apptainer/Singularity.
- `-profile conda`
  Системы, где проще управлять инструментами через Conda environments.

Важно: текущая архитектура не описывается как "Docker-контейнер, внутри которого Apptainer запускает `.sif`". Nextflow использует выбранный runtime напрямую:

- в `docker` профиле — Docker containers
- в `singularity` профиле — Singularity/Apptainer containers
- в `conda` профиле — `environment.yml` у модулей

Каталог `containers/` нужен для локальных кастомных модулей и сборки соответствующих runtime-образов, а не как обязательный слой вложенной контейнеризации.

## Конфигурация

### Основные параметры

| Параметр | Значение по умолчанию | Назначение |
|---|---|---|
| `--input` | `null` | CSV samplesheet для локальных FASTQ |
| `--sra_ids` | `null` | Текстовый файл со списком SRA accession |
| `--outdir` | `./results2` | Корневая директория результатов |
| `--reference` | `assets/h37rv.fa` | Референсный FASTA |
| `--gbk` | `assets/h37rv.gbk` | GenBank-файл H37Rv |
| `--snpeff_db` | `h37rv_custom` | Имя базы SnpEff |
| `--snpeff_data_dir` | `assets/SNPEFF_ANNOTATION/data` | Каталог данных SnpEff |
| `--multiqc_config` | `assets/multiqc/multiqc_config.yaml` | Конфиг MultiQC |
| `--mode` | `copy` | `publishDir` mode |
| `--min_align_pct` | `90` | Минимальный процент выравненных ридов |
| `--min_median` | `30` | Минимальное медианное покрытие |
| `--skip_qc` | `false` | Не запускать FastQC |
| `--skip_kraken` | `false` | Не запускать Kraken/Bracken |
| `--skip_multiqc` | `false` | Не строить MultiQC |
| `--skip_final_reports` | `false` | Не строить `Reports/general` и `Reports/tb-platform` |
| `--skip_snp_matrix` | `false` | Не строить cohort SNP matrix |

### Kraken2 / Bracken

Поддерживаются одна или две Kraken DB:

| Параметр | Назначение |
|---|---|
| `--kraken2_db` | Первая Kraken2 DB |
| `--kraken2_db_label` | Метка первой DB |
| `--kraken2_db_2` | Вторая Kraken2 DB |
| `--kraken2_db_label_2` | Метка второй DB |

Если label не задан, он выводится автоматически из имени каталога базы.

### Ресурсы по label

Определены в `conf/base.config`:

| Label | CPU | Memory | `maxForks` |
|---|---|---|---|
| `process_single` | 1 | 4 GB | 12 |
| `process_low` | 2 | 6 GB | 12 |
| `process_medium` | 6 | 8 GB | 6 |
| `process_high` | 6 | 8 GB | 2 |

## Примеры запуска

### FASTQ + Docker

```bash
python make_samplesheet.py -i data -o run.csv

nextflow run main.nf \
  -profile docker \
  --input run.csv \
  --outdir results \
  -resume
```

### FASTQ + Conda

```bash
nextflow run main.nf \
  -profile conda \
  --input run.csv \
  --outdir results \
  -resume
```

### SRA

```bash
printf "SRR32010433\nERR15166664\n" > sra_ids.txt

nextflow run main.nf \
  -profile docker \
  --sra_ids sra_ids.txt \
  --outdir results \
  -resume
```

### Kraken с одной базой

```bash
nextflow run main.nf \
  -profile docker \
  --input run.csv \
  --kraken2_db /path/to/kraken_db \
  --kraken2_db_label ALL \
  --outdir results \
  -resume
```

### Kraken с двумя базами

```bash
nextflow run main.nf \
  -profile docker \
  --input run.csv \
  --kraken2_db /path/to/db1 \
  --kraken2_db_label ALL \
  --kraken2_db_2 /path/to/db2 \
  --kraken2_db_label_2 ONLY_MYCOBACTERIUM \
  --outdir results \
  -resume
```

### Singularity / Apptainer

```bash
nextflow run main.nf \
  -profile singularity \
  --input run.csv \
  --outdir results \
  -resume
```

### Kubernetes

```bash
nextflow run main.nf \
  -profile docker \
  -c k8s.config \
  --input run.csv \
  --outdir results \
  -resume
```

## Batch-режим

Для длинных запусков по большим samplesheet используйте `run_batches.sh`.

Пример:

```bash
./run_batches.sh \
  --pipeline /home/zerg/git/tb-lite \
  --input /data/run.csv \
  --batch-size 500 \
  --profile conda \
  --outdir /data/results_batches
```

### Что делает batch-режим

1. Делит входной файл на батчи.
2. Каждый батч запускает `main.nf` с:
   - `--skip_final_reports`
   - `--skip_multiqc`
   - `--skip_snp_matrix`
3. Все sample-level outputs складываются в один общий `outdir`.
4. После последнего успешного батча запускается `batch_reports.nf`, который собирает:
   - один общий `Reports/`
   - один общий `multiqc/`

### Batch + Kraken

```bash
./run_batches.sh \
  --pipeline /home/zerg/git/tb-lite \
  --input /data/run.csv \
  --batch-size 500 \
  --profile conda \
  --with-kraken \
  --kraken2_db /data/kraken_db \
  --kraken2_db_label ALL \
  --kraken2_db_2 /data/myco_db \
  --kraken2_db_label_2 ONLY_MYCOBACTERIUM \
  --outdir /data/results_batches
```

По умолчанию `run_batches.sh` использует профиль `conda`.

## Примечания

- Для `SnpEff` пайплайн ожидает каталог данных в `assets/SNPEFF_ANNOTATION/data`.
- Для custom базы по умолчанию используется `--snpeff_db h37rv_custom`.
- `ISMapper` работает только с paired-end reads.
- Параметр `--samples` сохранён как устаревший алиас к `--input`.
- Параметр `--multiqc` сохранён как устаревший алиас к `--multiqc_config`.
