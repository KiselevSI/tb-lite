/*
 * call_variant.nf : call_variants → merge → annotate
 * IN : bam_good  – tuple val(sample_name), path(bam)
 *      reference – path(ref.fa)
 * OUT: vcf_annot – tuple val(sample_name), path(annot.vcf.gz)
 */

include { FREEBAYES } from '../../modules/nf-core/freebayes/main'
include { BCFTOOLS_INDEX } from '../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_NORM } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_VIEW } from '../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_STATS } from '../../modules/nf-core/bcftools/stats/main'
include { BCFTOOLS_ANNOTATE } from '../../modules/nf-core/bcftools/annotate/main'
include { SNPEFF_SNPEFF } from '../../modules/nf-core/snpeff/snpeff/main'

process COUNT_VARIANTS {
    tag "${meta.id}"
    label 'process_low'

    input:
        tuple val(meta), path(vcf), path(tbi)

    output:
        tuple val(meta), path("${meta.id}.var_count.txt")

    script:
        """
        gzip -cd ${vcf} | grep -cv '^#' > ${meta.id}.var_count.txt
        """
}

workflow CALLVAR {
    take:
        bam_good

    main:
        ref_meta = [id: file(params.reference).baseName]
        ref = Channel.value([ref_meta, file(params.reference)])
        ref_fai = Channel.value([ref_meta, file("${params.reference}.fai", checkIfExists: true)])
        snpeff_cache = Channel.value([[id: 'snpeff_cache'], file(params.snpeff_data_dir, checkIfExists: true)])
        chr_name = file(params.chr_name)

        call_variants_input = bam_good.map { sample_name, bam, bai ->
            [[id: sample_name], bam, bai, [], [], []]
        }

        FREEBAYES(
            call_variants_input,
            ref,
            ref_fai,
            [[], []],
            [[], []],
            [[], []]
        )

        BCFTOOLS_INDEX(FREEBAYES.out.vcf)

        freebayes_indexed = FREEBAYES.out.vcf
            .join(BCFTOOLS_INDEX.out.tbi)
            .map { meta, vcf_file, vcf_index -> [meta, vcf_file, vcf_index] }

        BCFTOOLS_NORM(freebayes_indexed, ref)

        norm_vcf = BCFTOOLS_NORM.out.vcf
            .join(BCFTOOLS_NORM.out.tbi)
            .map { meta, vcf_file, vcf_index -> [meta, vcf_file, vcf_index] }

        BCFTOOLS_VIEW(norm_vcf, [], [], [])

        filtered_vcf = BCFTOOLS_VIEW.out.vcf
            .join(BCFTOOLS_VIEW.out.tbi)
            .map { meta, vcf_file, vcf_index -> [meta, vcf_file, vcf_index] }

        BCFTOOLS_STATS(
            filtered_vcf,
            [[], []],
            [[], []],
            [[], []],
            [[], []],
            [[], []]
        )

        COUNT_VARIANTS(filtered_vcf)

        other_count = filtered_vcf
            .join(COUNT_VARIANTS.out)
            .map { meta, vcf_file, vcf_index, count_file ->
                tuple(meta.id, vcf_file, vcf_index, count_file)
            }

        vcf = filtered_vcf.map { meta, vcf_file, vcf_index ->
            tuple(meta.id, vcf_file, vcf_index)
        }

        if (params.snpeff_db == "Mycobacterium_tuberculosis_h37rv") {
            renamed_vcf_input = filtered_vcf.map { meta, vcf_file, vcf_index ->
                [meta, vcf_file, vcf_index, [], [], [], [], chr_name]
            }
            BCFTOOLS_ANNOTATE(renamed_vcf_input)
            snpeff_input = BCFTOOLS_ANNOTATE.out.vcf
        } else {
            snpeff_input = filtered_vcf.map { meta, vcf_file, vcf_index ->
                [meta, vcf_file]
            }
        }

        SNPEFF_SNPEFF(snpeff_input, params.snpeff_db, snpeff_cache)

        bcftools_stats = BCFTOOLS_STATS.out.stats.map { meta, stats_file -> stats_file }
        ann = SNPEFF_SNPEFF.out.vcf.map { meta, ann_vcf ->
            tuple(meta.id, ann_vcf)
        }

    emit:
        vcf
        other_count
        bcftools_stats
        ann
}
