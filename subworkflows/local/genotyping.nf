/*
 * genotyping.nf : spoligotyping + is6110 + rd
 * IN : bam_good  – tuple val(sample_name), path(bam)
 *      reference – path(ref.fa)
 * OUT: spol_table – path(<sample>.tsv)          — Sample/SpolBin/Spol8
 *      spol_log   – path(<sample>.log)          — лог SpoTyping (спейсеры, min/rmin)
 *      is6110     – tuple val(meta), path(results/*) — вывод ISMapper (только paired)
 *      del        – tuple val(sample_name), path(<sample>.rd.tsv)
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
        tblg_table = TBLG(vcf).lg

        paired_reads = trimmed_reads
            .filter { _id, files -> files.size() == 2 }

        PREPARE_ISMAPPER_READS(paired_reads)

        ismapper_reads = PREPARE_ISMAPPER_READS.out.reads
            .map { sample_name, files ->
                [[id: sample_name], files.sort { it.name }, file(params.gbk), file(params.is6110)]
            }

        is6110 = ISMAPPER(ismapper_reads).results

        SPOTYPING(trimmed_reads)
        spol_table = SPOTYPING.out.table
        spol_log = SPOTYPING.out.log

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
        spol_log
        is6110
        dr
        del

}
