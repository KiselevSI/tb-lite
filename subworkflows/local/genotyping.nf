/*
 * genotyping.nf : spoligotyping + is6110 + rd
 * IN : bam_good  – tuple val(sample_name), path(bam)
 *      reference – path(ref.fa)
 * OUT: spoligo   – tuple val(sample_name), path(spoligotyping.txt)
 *      is6110    – tuple val(sample_name), path(is6110.bed)
 *      rd        – tuple val(sample_name), path(read_depth.tab)
 */

include { ISMAPPER } from '../../modules/nf-core/ismapper/main'
include { SPOTYPING } from '../../modules/local/genotyping/spotyping/main'
include { MOSDEPTH } from '../../modules/nf-core/mosdepth/main'
include { RD } from '../../modules/local/genotyping/rd/main'
include { TBLG } from '../../modules/local/genotyping/tblg/main'
include { TB_PROFILER_DR } from '../../modules/local/genotyping/tb_profiler_dr/main'

process PREPARE_ISMAPPER_READS {
    tag "${sample_name}"
    label 'process_low'

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("*.fastq.gz"), emit: reads

    script:
    def files = reads instanceof List ? reads : [reads]
    """
    ln -sf ${files[0]} ${sample_name}_R1.fastq.gz
    ln -sf ${files[1]} ${sample_name}_R2.fastq.gz
    """
}

workflow GENOTYPE {
    take:
        trimmed_reads
        vcf
        bam_good

    main:
        tblg_table = TBLG(vcf)

        paired_reads = trimmed_reads
            .filter { _id, files -> files.size() == 2 }

        PREPARE_ISMAPPER_READS(paired_reads)

        ismapper_reads = PREPARE_ISMAPPER_READS.out.reads
            .map { sample_name, files ->
                [[id: sample_name], files.sort { it.name }, file(params.ref_gbk), file(params.is6110)]
            }

        ISMAPPER(ismapper_reads)

        spol_table = SPOTYPING(trimmed_reads).table

        ref_meta = [id: file(params.reference).baseName]
        ref = Channel.value([ref_meta, file(params.reference)])

        mosdepth_input = bam_good.map { sample_name, bam, bai ->
            [[id: sample_name], bam, bai, []]
        }

        MOSDEPTH(mosdepth_input, ref)

        bed = MOSDEPTH.out.per_base_bed.map { meta, bed_file ->
            tuple(meta.id, bed_file)
        }

        rd_db = Channel.value(file(params.rd_db))

        del = RD(bed, rd_db).rd

        TB_PROFILER_DR(bam_good)
        dr = TB_PROFILER_DR.out.json

    emit:
        tblg_table
        spol_table
        dr
        del

}
