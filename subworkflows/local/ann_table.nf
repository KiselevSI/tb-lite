
include { MAKE_TABLE } from '../../modules/local/ann_table/make_table/main'
include { BCFTOOLS_MERGE as ANN_BCFTOOLS_MERGE } from '../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_ANNOTATE as ANN_BCFTOOLS_ANNOTATE } from '../../modules/nf-core/bcftools/annotate/main'
include { SNPEFF_SNPEFF as ANN_SNPEFF } from '../../modules/nf-core/snpeff/snpeff/main'
include { BCFTOOLS_VIEW as ANN_BCFTOOLS_VIEW } from '../../modules/nf-core/bcftools/view/main'
include { POST_PROCESS_TABLE } from '../../modules/local/ann_table/post_process_table/main'

process FILTER_ANN_TABLE_INPUT {
    tag "${sample_name}"
    label 'process_low'

    input:
        tuple val(sample_name), path(vcf), path(tbi), path(count_file)

    output:
        tuple val(sample_name), path("${sample_name}.ann_input.vcf.gz"), path("${sample_name}.ann_input.vcf.gz.tbi"), optional: true, emit: kept

    script:
        """
        count_value=\$(tr -d '[:space:]' < ${count_file})

        if [ "\${count_value}" -lt 5000 ]; then
            ln -sf ${vcf} ${sample_name}.ann_input.vcf.gz
            ln -sf ${tbi} ${sample_name}.ann_input.vcf.gz.tbi
        fi
        """
}

workflow ANN_TABLE {
    take:
        vcfs
        ref_fai

    main:
        ref = ref_fai.map { meta, fai ->
            [meta, file(params.reference), fai]
        }
        snpeff_cache = Channel.value([[id: 'snpeff_cache'], file(params.snpeff_data_dir, checkIfExists: true)])

        FILTER_ANN_TABLE_INPUT(vcfs)

        filtered = FILTER_ANN_TABLE_INPUT.out.kept
            .collect(flat: false)
            .map { records ->
                if (records.size() <= 1) {
                    log.warn("Skipping SNP matrix generation: need at least 2 samples after ANN_TABLE filtering, got ${records.size()}.")
                }
                records
            }
            .filter { records -> records.size() > 1 }

        merge_input = filtered.map { records ->
            def cohort_vcfs = records.collect { it[1] }
            def cohort_tbis = records.collect { it[2] }
            tuple([id: 'cohort'], cohort_vcfs, cohort_tbis, [])
        }

        ANN_BCFTOOLS_MERGE(merge_input, ref)

        if (params.snpeff_db == "Mycobacterium_tuberculosis_h37rv") {
            chr_name = file(params.chr_name)
            ann_rename_input = ANN_BCFTOOLS_MERGE.out.vcf.map { meta, vcf ->
                [meta, vcf, [], [], [], [], [], chr_name]
            }
            ANN_BCFTOOLS_ANNOTATE(ann_rename_input)
            ann_snpeff_input = ANN_BCFTOOLS_ANNOTATE.out.vcf
        } else {
            ann_snpeff_input = ANN_BCFTOOLS_MERGE.out.vcf
        }

        ANN_SNPEFF(ann_snpeff_input, params.snpeff_db, snpeff_cache)

        ann_view_input = ANN_SNPEFF.out.vcf.map { meta, vcf ->
            [meta, vcf, []]
        }

        ANN_BCFTOOLS_VIEW(ann_view_input, [], [], [])

        table_input = ANN_BCFTOOLS_VIEW.out.vcf
            .join(ANN_BCFTOOLS_VIEW.out.tbi)
            .map { meta, vcf, tbi -> tuple(vcf, tbi) }

        table = MAKE_TABLE(table_input)
        feature_table = Channel.fromPath(params.h37rv_feature_table)
        final_table = POST_PROCESS_TABLE(table, feature_table)

    emit:
        clean = final_table.clean
}
