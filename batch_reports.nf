nextflow.enable.dsl=2

include { REPORTS } from './subworkflows/local/reports'

workflow BATCH_REPORTS {
    main:
        wgs_metrics = Channel.fromPath("${params.outdir}/stats/picard/wgs/*/*.CollectWgsMetrics.coverage_metrics", checkIfExists: true)
        align_metrics = Channel.fromPath("${params.outdir}/stats/picard/alignment/*/*.txt", checkIfExists: true)
        fastp_reports = Channel.fromPath("${params.outdir}/fastp/*/*.fastp.json", checkIfExists: true)
        fastqc_reports = Channel.fromPath("${params.outdir}/fastqc/*/*_fastqc.zip", checkIfExists: true)
        bcftools_stats = Channel.fromPath("${params.outdir}/stats/bcftools/*.bcftools_stats.txt", checkIfExists: true)
        samtools_stats = Channel.fromPath("${params.outdir}/stats/samtools/stats/*/*.stats", checkIfExists: true)
        samtools_flagstat = Channel.fromPath("${params.outdir}/stats/samtools/flagstat/*/*.flagstat", checkIfExists: true)
        tbmix = Channel.fromPath("${params.outdir}/tb-mix/*.mix.tsv", checkIfExists: true)
        spotyping = Channel.fromPath("${params.outdir}/spotyping/*/*.tsv", checkIfExists: true)
        tblg_table = Channel.fromPath("${params.outdir}/lineage/*.lg.tsv", checkIfExists: true)
        drugs = Channel.fromPath("${params.outdir}/tb-profiler/drug-resist/*/results/*.results.json", checkIfExists: true)
        del = Channel
            .fromPath("${params.outdir}/rd/*/*.novel_rd.tsv", checkIfExists: true)
            .map { path -> tuple(path.name.replaceFirst(/\.novel_rd\.tsv$/, ''), path) }

        kraken_reports = (!params.skip_kraken && params.kraken2_db)
            ? Channel.fromPath("${params.outdir}/kraken2/kraken2/*/*/*.kraken2.report.txt", checkIfExists: true)
            : Channel.empty()

        kraken_combined = (!params.skip_kraken && params.kraken2_db)
            ? Channel
                .fromPath("${params.outdir}/kraken2/combined/*.all_samples.txt", checkIfExists: true)
                .map { path -> tuple([id: path.baseName], path) }
            : Channel.empty()

        REPORTS(
            wgs_metrics,
            align_metrics,
            fastp_reports,
            fastqc_reports,
            bcftools_stats,
            samtools_stats,
            samtools_flagstat,
            kraken_reports,
            kraken_combined,
            tbmix,
            spotyping,
            tblg_table,
            drugs,
            del
        )
}

workflow {
    BATCH_REPORTS()
}
