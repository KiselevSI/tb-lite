/*
 * trimming.nf : fastp
 * IN : tuple val(sample_name), path(reads)
 * OUT: trimmed_reads – tuple val(sample_name), path(*.trimmed.fq.gz)
 */
include { FASTP } from '../../modules/local/trimming/fastp/main'
include { REPAIR_READS } from '../../modules/local/trimming/repair_reads/main'
workflow TRIMMING {
    take:
        raw_reads                                     // (sample_name, reads)

    main:
        //repair_reads = REPAIR_READS(raw_reads).repaired_reads
        trimmed_reads = FASTP(raw_reads).trimmed_reads

    emit:
        trimmed_reads
}
