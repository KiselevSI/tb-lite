/*
 * trimming.nf : fastp
 * IN : tuple val(sample_name), path(reads)
 * OUT: trimmed_reads – tuple val(sample_name), path(*.trimmed.fq.gz)
 */
include { FASTP } from '../modules/trimming/fastp'
include { REPAIR_READS } from '../modules/trimming/repair_reads'
workflow TRIMMING {
    take:
        raw_reads                                     // (sample_name, reads)

    main:
        //repair_reads = REPAIR_READS(raw_reads).repaired_reads
        trimmed_reads = FASTP(raw_reads).trimmed_reads

    emit:
        trimmed_reads
}
