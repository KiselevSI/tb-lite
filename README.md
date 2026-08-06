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

| Файл | Содержимое | Импорт в TB Platform |
| --- | --- | --- |
| `general.tsv` | Метрики покрытия и выравнивания | `scripts/import_new_core_data.py` |
| `filter.tbmix.tsv` | TB-Mix + линии после фильтрации | — |
| `drug_resist.xlsx` | Лекарственная устойчивость | `backend/scripts/import_drug_resist_xlsx.py` |
| `drug_resist_and_uncertain.xlsx` | То же плюс uncertain-варианты | — |
| `spotyping.total.tsv` | `Sample`/`SpolBin`/`Spol8` | `general_spoligo` |
| `spotyping.full.tsv` | То же плюс `MinReads`/`RminReads` и SIT/клада/география из SpolDB4 | справочно |
| `spoligo_spacer_counts.tsv` | Число ридов на каждый из 43 спейсеров | `general_spoligo_spacers` |
| `rd.tsv` | Known + novel RD-делеции, 19 колонок (`rd_scan.py`) | `scripts/import_deletions_full.py` |

Процесс `IS6110_TABLES` публикует в `Reports/tb-platform/is6110/` нормализованные
таблицы вставок IS6110 — каталог целиком принимает
`backend/scripts/import_is6110_tsv.py --input-dir`:

- `is_element.tsv`, `is6110_site.tsv`, `ismapper_run.tsv`,
  `sample_is6110_site.tsv`, `is6110_removed_hit.tsv`, `is6110_sample_summary.tsv`
- диагностика: `import_stats.tsv`, `missing_samples.tsv`, `duplicate_runs.tsv`,
  `is6110_import_warnings.tsv`

ISMapper запускается только по paired-образцам, поэтому для набора целиком из
single-end данных таблицы публикуются пустыми (только шапки).

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

По умолчанию sample ID читается из VCF header. Для больших наборов это надежно, но медленно, потому что надо открыть каждый VCF. Если sample ID можно брать из имени файла, используйте быстрый режим:

```bash
python make_snp_matrix_csv.py \
  -i ../VCF/ /data6/bio/MolGenMicro/TBGenoPipe/results/VCF2/VCF \
  -o vcf_samples.csv \
  --sample-source filename
```

Для `--sample-source filename` параметр `--jobs` распараллеливает обход входных директорий и запись временных shard-файлов. Прогресс поиска и обработки печатается в stderr:

```bash
python make_snp_matrix_csv.py \
  -i ../VCF/ /data6/bio/MolGenMicro/TBGenoPipe/results/VCF2/VCF \
  -o vcf_samples.csv \
  --sample-source filename \
  --jobs 16
```

Если нужен именно sample ID из header, включите параллельное чтение:

```bash
python make_snp_matrix_csv.py \
  -i ../VCF/ /data6/bio/MolGenMicro/TBGenoPipe/results/VCF2/VCF \
  -o vcf_samples.csv \
  --sample-source header \
  --jobs 16
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
   - один общий `Reports/` (включая `Reports/tb-platform/snp/`)
   - один общий `multiqc/`

### Только SNP-таблицы для готового прогона

Если каталог с результатами уже есть, а SNP-таблиц в нём нет, их можно
собрать отдельно — без остальных отчётов:

```bash
nextflow run batch_reports.nf -entry SNP_DB_ONLY \
  -profile docker \
  --outdir /data/results_batches
```

Если `snp_db/shards/` пуст, шарды соберутся из `annotate_vcf/` пачками по 500.
Отдельная точка входа нужна потому, что полный `batch_reports.nf` объявляет
входные каналы с `checkIfExists: true` и падает целиком, если в каталоге
прогона пуст хотя бы один из его глобов.

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

## SNP-таблицы для TB Platform

Пайплайн генерирует файлы, которыми образцы добавляются на сайт. Итог —
в `<outdir>/Reports/tb-platform/snp/`:

| Файл | Назначение |
|---|---|
| `snp_sites.tsv.gz` | каталог аннотированных сайтов; `site_id` выдаёт БД по `site_key` |
| `sample_snp_alleles.tsv.gz` | аллели образцов — SNP matrix TSV и поиск по SNP |
| `vcf_table.tsv.gz` | плоская таблица `sample_id/pos/alt` — поиск и clade-признаки |
| `import_snp.sql` | psql-скрипт загрузки |
| `snp_db_manifest.tsv` | счётчики строк для сверки |
| `samples.txt` | список загружаемых образцов |

Источник — per-sample `annotate_vcf/<sample>/<sample>.annotated.ann.vcf` после snpEff.
Пишется **всё содержимое VCF**: SNP, инделы, complex и mnp, без фильтров и без масок.

`site_key` — sha1 от `chrom, pos, ref, alt` и полей ANN. Формула совпадает
с продовой таблицей `snp_sites`, поэтому новые данные схлопываются
с существующими через `ON CONFLICT (site_key)`.

Отключается флагом `--skip_snp_db`.

### Как это работает в batch-режиме

1. Каждый батч `main.nf` пишет свои шарды в `<outdir>/snp_db/shards/`
   (одна задача на батч, тег из `--batch_tag`, который `run_batches.sh`
   уже передаёт). Перезапуск батча перезаписывает те же файлы, дублей не будет.
2. Финальный `batch_reports.nf` склеивает все шарды в три таблицы,
   дедуплицируя `snp_sites` по `site_key` с максимальным `QUAL`.

`run_batches.sh` менять не нужно: `--batch_tag batch_${i}` он уже передаёт,
а финальный `batch_reports.nf` подхватывает готовые шарды.

Не запускайте batch-режим с `--mode link` или `--mode symlink`: `run_batches.sh`
чистит `work/` после каждого батча, и шарды-симлинки станут битыми. По умолчанию
`mode = copy`, этого и держитесь.

Если прогон уже отработал без шардов, их можно собрать задним числом
из `annotate_vcf/` — см. `-entry SNP_DB_ONLY` ниже.

### Импорт в базу

Порядок важен: `general.tsv` должен быть импортирован раньше — на
`sample_snp_alleles.sample_id` висит внешний ключ на `general."ID"`.

```bash
psql "$DATABASE_URL" -f <outdir>/Reports/tb-platform/snp/import_snp.sql
```

`\copy` выполняется на стороне клиента, поэтому каталог с `.tsv.gz` должен быть
виден тому процессу, который запускает psql. Если psql запускается в контейнере,
каталог надо смонтировать и перегенерировать SQL с путём внутри контейнера:

```bash
python bin/write_import_snp_sql.py --data-dir /mnt/snp -o import_snp.sql
```

Повторный запуск импорта безопасен: `snp_sites` схлопывается по `site_key`,
`sample_snp_alleles` — по первичному ключу, `vcf_table` перезаписывается по
списку входящих образцов, профили обновляются.

Последним шагом тот же SQL достраивает `sample_snp_profiles` по маскам из
`snp_masks` / `snp_mask_regions`. Если таблица масок пуста, шаг пропускается
с предупреждением — маски заливаются отдельно
(`tb-platform/deploy/rebuild_snp_profiles_masked.sql`).

## Примечания

- Для `SnpEff` пайплайн ожидает каталог данных в `assets/SNPEFF_ANNOTATION/data`.
- Для custom базы по умолчанию используется `--snpeff_db h37rv_custom`.
- `ISMapper` работает только с paired-end reads.
- Параметр `--samples` сохранён как устаревший алиас к `--input`.
- Параметр `--multiqc` сохранён как устаревший алиас к `--multiqc_config`.
