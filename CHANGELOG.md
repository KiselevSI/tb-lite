# Changelog

## Unreleased

### Changed

- **BREAKING.** RD-детекция переведена с `bin/rd.py` на `bin/rd_scan.py`. Процесс
  `RD` теперь публикует один `<sample>.rd.tsv` (19 колонок, known + novel) вместо
  пары `novel_rd.tsv` / `known_rd.tsv`. Соответственно `Reports/tb-platform/rd.tsv`
  из 9-колоночного CSV стал 19-колоночным TSV — формат, который принимает
  `scripts/import_deletions_full.py` в TB Platform. Образ обновлён до
  `tb-lite/rd-scanner:2.0`, зависимости numpy/pandas больше не нужны.
- `bin/deletions_to_csv.py` удалён: per-sample таблицы склеиваются
  `bin/concat_tables.py --keep-header`.

### Added

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
