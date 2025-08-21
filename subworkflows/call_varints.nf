/*
 * call_variant.nf : call_variants → merge → annotate
 * IN : bam_good  – tuple val(sample_name), path(bam)
 *      reference – path(ref.fa)
 * OUT: vcf_annot – tuple val(sample_name), path(annot.vcf.gz)
 */

include { BCFTOOLS_CALL_VARIANTS } from '../modules/call_variants/bcftools'
include { BCFTOOLS_STATS } from '../modules/call_variants/bcftools_stats'
include { RENAME_CHR } from '../modules/call_variants/rename_chr'
include { SNPEFF_ANNOTATE_VCF } from '../modules/call_variants/snpeff'

workflow CALLVAR {
    take:
        bam_good

    main:
        chr_name = Channel.value(file(params.chr_name))
        ref = Channel.value(file(params.reference))
        call_variants = BCFTOOLS_CALL_VARIANTS(bam_good, ref)
        vcf = call_variants.other
        bcftools_stats = BCFTOOLS_STATS(vcf)
        vcfs = call_variants.vcfs
        vcf_renamed = RENAME_CHR(vcf, chr_name)
        vcf_annotated = SNPEFF_ANNOTATE_VCF(vcf_renamed).ann

    emit:
        vcf
        vcfs
        bcftools_stats
        vcf_annotated
}
