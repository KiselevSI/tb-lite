# Changelog

## Unreleased

### Changed

- Таблицы `Reports/tb-platform/general.tsv`, `tbmix.total.tsv` и
  `filter.tbmix.tsv` теперь содержат только образцы, дошедшие до вызова
  вариантов. Раньше метрики покрытия и TB-Mix считались до фильтра качества,
  и отбракованные образцы (нулевое покрытие) попадали в таблицы для базы.
  Фильтрация — по списку `stats/bcftools/*.bcftools_stats.txt`, скрипт
  `bin/filter_table_by_samples.py`. QC-отчёты `Reports/general/` по-прежнему
  показывают все образцы.
- **BREAKING.** RD-детекция переведена с `bin/rd.py` на `bin/rd_scan.py`. Процесс
  `RD` теперь публикует один `<sample>.rd.tsv` (19 колонок, known + novel) вместо
  пары `novel_rd.tsv` / `known_rd.tsv`. Соответственно `Reports/tb-platform/rd.tsv`
  из 9-колоночного CSV стал 19-колоночным TSV — формат, который принимает
  `scripts/import_deletions_full.py` в TB Platform. Образ обновлён до
  `tb-lite/rd-scanner:2.0`, зависимости numpy/pandas больше не нужны.
- `bin/deletions_to_csv.py` удалён: per-sample таблицы склеиваются
  `bin/concat_tables.py --keep-header`.

### Added

- `Reports/tb-platform/snp/` — таблицы VCF для TB Platform: `snp_sites.tsv.gz`,
  `sample_snp_alleles.tsv.gz`, `vcf_table.tsv.gz` плюс `import_snp.sql`,
  `snp_db_manifest.tsv` и `samples.txt`. Новые процессы `SNP_DB_SHARDS`
  (шард на батч) и `SNP_DB_TABLES` (склейка), скрипты
  `bin/vcf_to_snp_shards.py`, `bin/merge_snp_shards.sh`,
  `bin/write_import_snp_sql.py`. Отключается флагом `--skip_snp_db`.
  Точка входа `batch_reports.nf -entry SNP_DB_ONLY` собирает только эти
  таблицы, при необходимости восстанавливая шарды из `annotate_vcf/`.
- `Reports/tb-platform/is6110/` — нормализованные таблицы вставок IS6110
  (новый процесс `IS6110_TABLES`, скрипт `bin/build_is6110_tables.py`).
  Раньше вывод ISMapper оставался только в per-sample каталогах.
- `Reports/tb-platform/spoligo_spacer_counts.tsv` и `spotyping.full.tsv` —
  число ридов на 43 спейсера, пороги `min`/`rmin` и SIT/клада/география из
  SpolDB4 (`bin/build_spoligo_table.py`, справочник
  `assets/spoldb4/spoldb4_reference.tsv`, параметр `--spoldb4`).
  `spotyping.total.tsv` формат не меняет.

## v1.0.0 - 2026-03-30

### Added

- Initial nf-core-compatible release
- 8-stage pipeline: trimming, QC, mapping, filtering, variant calling, genotyping, annotation, reports
- Support for paired-end and single-end reads
- Three variant callers: FreeBayes (default), GATK, bcftools
- Drug resistance profiling via TB-Profiler
- Spoligotyping, lineage classification, IS6110 mapping, region of difference analysis
- Docker and Singularity container support
- Local and Kubernetes execution profiles
- Batch execution script for large-scale runs (run_batches.sh)
