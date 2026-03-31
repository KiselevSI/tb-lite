# Changelog

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
