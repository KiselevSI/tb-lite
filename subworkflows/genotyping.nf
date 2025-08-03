/*
 * genotyping.nf : spoligotyping + is6110 + rd
 * IN : bam_good  – tuple val(sample_name), path(bam)
 *      reference – path(ref.fa)
 * OUT: spoligo   – tuple val(sample_name), path(spoligotyping.txt)
 *      is6110    – tuple val(sample_name), path(is6110.bed)
 *      rd        – tuple val(sample_name), path(read_depth.tab)
 */

include { ISMAPPER } from '../modules/genotyping/ismapper'
include { SPOTYPING } from '../modules/genotyping/spotyping'
include { MOSDEPTH } from '../modules/genotyping/mosdepth'
include { RD } from '../modules/genotyping/rd'
include { TBLG } from '../modules/genotyping/tblg'
include { TB_PROFILER_DR } from '../modules/genotyping/tb_profiler'

workflow GENOTYPE {
    take:
        paired_reads
        trimmed_reads
        vcf
        bam_good

    main:
        tblg_table = TBLG(vcf)


        is6110 = Channel.value(file(params.is6110))
        ref_gbk = Channel.value(file(params.ref_gbk))


        ISMAPPER(paired_reads, is6110, ref_gbk)
        spol_table = SPOTYPING(trimmed_reads).table

        bed = MOSDEPTH(bam_good).bed

        rd_db = Channel.value(file(params.rd_db))

        RD(bed, rd_db)

        dr = TB_PROFILER_DR(vcf).results

    emit:
        tblg_table
        spol_table
        dr

}
