include { MULTIQC } from '../../modules/nf-core/multiqc/main'
include { FINAL_TABLE } from '../../modules/local/reports/final_table/main'
include { TB_PLATFORM_TABLES } from '../../modules/local/reports/tb_platform_tables/main'
include { IS6110_TABLES } from '../../modules/local/reports/is6110_tables/main'

workflow REPORTS {
    take:
        wgs_metrics
        align_metrics
        fastp_reports
        fastqc_reports
        bcftools_stats
        samtools_stats
        samtools_flagstat
        kraken_reports
        kraken_combined
        tbmix
        spotyping
        spotyping_logs
        tblg_table
        drugs
        del
        is6110

    main:
        skip_multiqc = params.skip_multiqc || params.skip_reports
        skip_final_reports = params.skip_final_reports || params.skip_reports
        final_table = Channel.empty()
        versions = Channel.empty()
        rd = del.map { _sample_name, path -> path }  // извлекаем только пути из канала del

        // ISMAPPER эмитит path("results/*") — это каталог results/<sample>, а не
        // отдельные файлы. Разворачиваем его и берём только нужные таблицы:
        // <sample>__<ref>_table.txt (вызовы) и <sample>__<ref>_table (удалённые хиты).
        is6110_tables = is6110
            .map { _meta, entries -> entries }
            .flatten()
            .flatMap { entry ->
                entry.isDirectory() ? files("${entry}/**/*__*_table*") : [entry]
            }
            .filter { it.name ==~ /.+__.+_table(\.txt)?$/ }
        kraken_combined_tables = (!params.skip_kraken && params.kraken2_db)
            ? kraken_combined.map { meta, path -> path }.collect()
            : Channel.value([])

        gbk = Channel.value(file(params.gbk))
        spoldb4 = Channel.value(file(params.spoldb4))

        if (!skip_multiqc) {
            multiqc_files = wgs_metrics.mix(align_metrics)
                .mix(fastp_reports)
                .mix(fastqc_reports)
                .mix(bcftools_stats)
                .mix(samtools_stats)
                .mix(samtools_flagstat)
                .mix(kraken_reports)

            cfg_path = params.multiqc_config ?: params.multiqc
            cfg = cfg_path ? file(cfg_path) : []

            MULTIQC(
                multiqc_files.collect().map { files ->
                    [[id: 'multiqc'], files, cfg, [], [], []]
                }
            )
            versions = MULTIQC.out.versions.map { _process_name, program, version ->
                tuple(program, version)
            }
        }

        if (!skip_final_reports) {
            TB_PLATFORM_TABLES(
                wgs_metrics.collect(),
                bcftools_stats.collect(),
                samtools_flagstat.collect(),
                tbmix.collect(),
                spotyping.collect(),
                spotyping_logs.collect(),
                tblg_table.collect(),
                drugs.collect(),
                rd.collect(),
                gbk,
                spoldb4
            )

            // ISMapper запускается только по paired-образцам — прогон целиком из
            // single-end данных не должен ронять отчёты.
            IS6110_TABLES(is6110_tables.collect().ifEmpty([]))

            final_table = FINAL_TABLE(
                wgs_metrics.collect(),
                fastp_reports.collect(),
                bcftools_stats.collect(),
                samtools_flagstat.collect(),
                tbmix.collect(),
                spotyping.collect(),
                tblg_table.collect(),
                drugs.collect(),
                gbk,
                kraken_combined_tables
            )
        }


    emit:
        final_table
        versions
}
