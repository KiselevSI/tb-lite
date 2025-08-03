/*
 * mapping.nf : bwa_index → run_mapping
 * IN : trimmed_reads – tuple val(sample_name), path(reads)
 *      reference     – path(ref.fa)
 * OUT: bam           – tuple val(sample_name), path(sorted.bam)
 */
include { BWA_INDEX } from '../modules/mapping/bwa_index'
include { BWA_PICARD } from '../modules/mapping/bwa_picard'
workflow MAPPING {
    take:
        trimmed_reads

    main:
        ref = Channel.value(file(params.reference))
        index = BWA_INDEX(ref).index                           // value-канал с индексом
        bam = BWA_PICARD(trimmed_reads, index, ref).bam

    emit:
        bam
}
