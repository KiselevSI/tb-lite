include { KRAKEN2_KRAKEN2 as KRAKEN2_DB1 }               from '../../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_DB2 }               from '../../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN as BRACKEN_DB1 }               from '../../modules/nf-core/bracken/bracken/main'
include { BRACKEN_BRACKEN as BRACKEN_DB2 }               from '../../modules/nf-core/bracken/bracken/main'
include { BRACKEN_COMBINEBRACKENOUTPUTS as COMBINE_DB1 } from '../../modules/nf-core/bracken/combinebrackenoutputs/main'
include { BRACKEN_COMBINEBRACKENOUTPUTS as COMBINE_DB2 } from '../../modules/nf-core/bracken/combinebrackenoutputs/main'
include { BRACKEN_ADD_UNCLASSIFIED as ADD_DB1 }          from '../../modules/local/kraken/bracken_add_unclassified'
include { BRACKEN_ADD_UNCLASSIFIED as ADD_DB2 }          from '../../modules/local/kraken/bracken_add_unclassified'

def sanitizeLabel(value) {
    value.toString().replaceAll(/[^A-Za-z0-9._-]+/, '_').replaceAll(/^_+|_+$/, '')
}

def deriveLabel(dbPath, explicitLabel) {
    def rawLabel = explicitLabel ?: file(dbPath).getBaseName()
    def label = sanitizeLabel(rawLabel)
    if (!label) {
        throw new IllegalArgumentException("Could not derive database label for '${dbPath}'")
    }
    label
}

workflow KRAKEN {
    take:
    trimmed_reads

    main:
    db1_label = deriveLabel(params.kraken2_db, params.kraken2_db_label)
    db2_label = params.kraken2_db_2 ? deriveLabel(params.kraken2_db_2, params.kraken2_db_label_2) : null

    if (db2_label && db1_label == db2_label) {
        throw new IllegalArgumentException("Database labels must be unique. Got '${db1_label}' for both databases.")
    }

    report = Channel.empty()
    bracken = Channel.empty()
    combined = Channel.empty()
    multiqc = Channel.empty()

    db1_reads = trimmed_reads.map { meta, fastqs ->
        tuple(meta + [db_label: db1_label], fastqs)
    }
    db1_path = Channel.value(file(params.kraken2_db, checkIfExists: true))

    KRAKEN2_DB1(db1_reads, db1_path, false, false)
    BRACKEN_DB1(KRAKEN2_DB1.out.report, db1_path)

    db1_bracken_kraken = BRACKEN_DB1.out.reports.join(KRAKEN2_DB1.out.report)
    ADD_DB1(db1_bracken_kraken)

    db1_all_bracken = ADD_DB1.out.tsv
        .map { meta, tsv -> tsv }
        .collect()

    COMBINE_DB1(
        db1_all_bracken.map { files -> [[id: 'all_samples', db_label: db1_label], files] }
    )

    report = report.mix(KRAKEN2_DB1.out.report)
    bracken = bracken.mix(ADD_DB1.out.tsv)
    combined = combined.mix(COMBINE_DB1.out.txt)
    multiqc = multiqc.mix(KRAKEN2_DB1.out.report.map { meta, kraken_report -> kraken_report })

    if (params.kraken2_db_2) {
        db2_reads = trimmed_reads.map { meta, fastqs ->
            tuple(meta + [db_label: db2_label], fastqs)
        }
        db2_path = Channel.value(file(params.kraken2_db_2, checkIfExists: true))

        KRAKEN2_DB2(db2_reads, db2_path, false, false)
        BRACKEN_DB2(KRAKEN2_DB2.out.report, db2_path)

        db2_bracken_kraken = BRACKEN_DB2.out.reports.join(KRAKEN2_DB2.out.report)
        ADD_DB2(db2_bracken_kraken)

        db2_all_bracken = ADD_DB2.out.tsv
            .map { meta, tsv -> tsv }
            .collect()

        COMBINE_DB2(
            db2_all_bracken.map { files -> [[id: 'all_samples', db_label: db2_label], files] }
        )

        report = report.mix(KRAKEN2_DB2.out.report)
        bracken = bracken.mix(ADD_DB2.out.tsv)
        combined = combined.mix(COMBINE_DB2.out.txt)
        multiqc = multiqc.mix(KRAKEN2_DB2.out.report.map { meta, kraken_report -> kraken_report })
    }

    emit:
    report
    bracken
    combined
    multiqc
}
