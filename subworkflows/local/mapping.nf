/*
 * mapping.nf : nf-core bwa + picard
 * IN : trimmed_reads – tuple val(sample_name), path(reads)
 * OUT: bam           – tuple val(sample_name), path(dedup.bam), path(dedup.bam.bai)
 */
include { BWA_INDEX }               from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM }                 from '../../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_FAIDX }          from '../../modules/nf-core/samtools/faidx/main'
include { PICARD_MARKDUPLICATES }   from '../../modules/nf-core/picard/markduplicates/main'

workflow MAPPING {
    take:
        trimmed_reads

    main:
        ref_meta = [id: file(params.reference).baseName]
        ref_ch = Channel.value([ref_meta, file(params.reference)])

        BWA_INDEX(ref_ch)

        ch_faidx_input = ref_ch.map { meta, fasta -> [meta, fasta, []] }
        SAMTOOLS_FAIDX(ch_faidx_input, false)

        ref_fai = SAMTOOLS_FAIDX.out.fai

        ch_reads_meta = trimmed_reads.map { sample_name, reads ->
            [[id: sample_name], reads]
        }

        BWA_MEM(
            ch_reads_meta,
            BWA_INDEX.out.index,
            ref_ch,
            true
        )

        ch_ref_fasta_fai = ref_fai.map { meta, fai ->
            [meta, file(params.reference), fai]
        }

        PICARD_MARKDUPLICATES(
            BWA_MEM.out.bam,
            ch_ref_fasta_fai
        )

        bam = PICARD_MARKDUPLICATES.out.bam
            .join(PICARD_MARKDUPLICATES.out.bai)
            .map { meta, bam_file, bai_file -> tuple(meta.id, bam_file, bai_file) }

    emit:
        bam
}
