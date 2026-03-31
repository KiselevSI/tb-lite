/*
 * qc.nf : fastqc
 * IN : tuple val(sample_name), path(reads)
 * OUT: fastqc_reports – tuple val(sample_name), path(*_fastqc.{html,zip})
 */
include { FASTQC } from '../../modules/local/qc/fastqc/main'

workflow QC {
    take:
        trimmed_reads

    main:
        reports = FASTQC(trimmed_reads)

    emit:
        reports
}
