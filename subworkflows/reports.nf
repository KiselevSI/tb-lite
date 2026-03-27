/*
 * reports.nf : сбор всех отчётов в одну папку/архив
 * IN : fastqc_reports – tuple val(sample_name), path(html/zip)
 *      tbmix          – tuple val(sample_name), path(tb_mix.txt)
 *      variant_stats  – tuple val(sample_name), path(bcftools_stats.txt)
 * OUT: multiqc_zip    – path(multiqc_report.zip)
 */
include { MULTIQC } from '../modules/reports/multiqc'   // условный пример
include { FINAL_TABLE } from '../modules/reports/final_table'   // условный пример
include { TB_PLATFORM_TABLES } from '../modules/reports/tb_platform_tables'   // условный пример

workflow REPORTS {
    take:
        wgs_metrics
        align_metrics
        fastqc_reports
        bcftools_stats
        samtools_stats
        samtools_flagstat
        tbmix
        spotyping
        tblg_table
        drugs
        del

    main:
        multiqc_files = wgs_metrics.mix(align_metrics)
            .mix(fastqc_reports)
            .mix(bcftools_stats)
            .mix(samtools_stats)
            .mix(samtools_flagstat)

        cfg = Channel.fromPath(params.multiqc)

        multiqc = MULTIQC(multiqc_files.collect(), cfg)

        rd = del.map { _sample_name, path -> path }  // извлекаем только пути из канала del

        gbk = Channel.value(file(params.gbk))

        TB_PLATFORM_TABLES(
            multiqc.report,
            tbmix.collect(),
            spotyping.collect(),
            tblg_table.collect(),
            drugs.collect(),
            rd.collect(),
            gbk
        )
        final_table = FINAL_TABLE(multiqc.report, tbmix.collect(), spotyping.collect(), tblg_table.collect(), drugs.collect())


    emit:
        final_table
}
