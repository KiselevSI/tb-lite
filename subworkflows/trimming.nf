/*
 * trimming.nf : fastp
 * IN : tuple val(sample_name), path(reads)
 * OUT: trimmed_reads – tuple val(sample_name), path(*.trimmed.fq.gz)
 */
include { FASTP } from '../modules/trimming/fastp'

workflow TRIMMING {
    take:
        raw_reads                                     // (sample_name, reads)

    main:
        trimmed_reads = FASTP(raw_reads).trimmed_reads

    emit:
        trimmed_reads
}
