

include { RENAME_CHR } from '../modules/vcf2table/rename_chr'
include { ANNOTATION_VCFS } from '../modules/vcf2table/snpeff'
include { MERGE_VCF } from '../modules/vcf2table/merge_vcf'
include { ANNOTATION_TABLE } from '../modules/vcf2table/annotation_table'

workflow VCF2TABLE {
    take:
        vcfs

    main:
        chr_name = Channel.value(file(params.chr_name))
        merged_vcf = MERGE_VCF(vcfs.collect())
        renamed_vcf = RENAME_CHR(merged_vcf, chr_name)
        merged_vcf_annotated = ANNOTATION_VCFS(renamed_vcf)

        feature_table = Channel.value(file(params.feature_table))
        ann_table = ANNOTATION_TABLE(vcfs.collect(), merged_vcf_annotated, feature_table)

    emit:
        ann_table
}