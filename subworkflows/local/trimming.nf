/*
 * trimming.nf : nf-core/fastp
 * IN : tuple val(meta), path(reads)
 * OUT: trimmed_reads – tuple val(meta), path(*.fastp.fastq.gz)
 */
include { FASTP } from '../../modules/nf-core/fastp/main'

workflow TRIMMING {
    take:
        raw_reads

    main:
        ch_fastp_input = raw_reads.map { meta, reads -> [meta, reads, []] }

        FASTP(ch_fastp_input, false, false, false)

        trimmed_reads = FASTP.out.reads
        trimmed_reads_legacy = FASTP.out.reads.map { meta, reads -> tuple(meta.id, reads) }
        fastp_json = FASTP.out.json.map { meta, json -> json }

    emit:
        trimmed_reads
        trimmed_reads_legacy
        fastp_json
}
