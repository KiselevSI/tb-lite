include { KRAKEN2_KRAKEN2 }                from '../../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN }                from '../../modules/nf-core/bracken/bracken/main'
include { BRACKEN_COMBINEBRACKENOUTPUTS }  from '../../modules/nf-core/bracken/combinebrackenoutputs/main'
include { BRACKEN_ADD_UNCLASSIFIED }       from '../../modules/local/kraken/bracken_add_unclassified'

workflow KRAKEN {
    take:
    trimmed_reads

    main:
    ch_kraken_db = Channel.value(file(params.kraken2_db, checkIfExists: true))

    KRAKEN2_KRAKEN2(trimmed_reads, ch_kraken_db, false, false)
    BRACKEN_BRACKEN(KRAKEN2_KRAKEN2.out.report, ch_kraken_db)

    ch_bracken_kraken = BRACKEN_BRACKEN.out.reports.join(KRAKEN2_KRAKEN2.out.report)
    BRACKEN_ADD_UNCLASSIFIED(ch_bracken_kraken)

    ch_all_bracken = BRACKEN_ADD_UNCLASSIFIED.out.tsv
        .map { meta, tsv -> tsv }
        .collect()

    BRACKEN_COMBINEBRACKENOUTPUTS(
        ch_all_bracken.map { files -> [[id: 'all_samples'], files] }
    )

    multiqc = KRAKEN2_KRAKEN2.out.report.map { meta, report -> report }

    emit:
    report = KRAKEN2_KRAKEN2.out.report
    bracken = BRACKEN_ADD_UNCLASSIFIED.out.tsv
    combined = BRACKEN_COMBINEBRACKENOUTPUTS.out.txt
    multiqc
}
