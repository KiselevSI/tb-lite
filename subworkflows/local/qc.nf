/*
 * qc.nf : nf-core/fastqc
 * IN : tuple val(meta), path(reads)
 * OUT: reports – path(*_fastqc.zip)
 */
include { FASTQC } from '../../modules/nf-core/fastqc/main'

workflow QC {
    take:
        trimmed_reads

    main:
        FASTQC(trimmed_reads)
        reports = FASTQC.out.zip.map { meta, zip -> zip }

    emit:
        reports
}
