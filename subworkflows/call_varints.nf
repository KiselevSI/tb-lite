/*
 * call_variant.nf : call_variants → merge → annotate
 * IN : bam_good  – tuple val(sample_name), path(bam)
 *      reference – path(ref.fa)
 * OUT: vcf_annot – tuple val(sample_name), path(annot.vcf.gz)
 */

include { BCFTOOLS_CALL_VARIANTS } from '../modules/call_variants/bcftools'
include { FREEBAYES_CALL_VARIANTS } from '../modules/call_variants/freebayes'
include { GATK_CALL_VARIANTS } from '../modules/call_variants/gatk'
include { BCFTOOLS_STATS } from '../modules/call_variants/bcftools_stats'
include { RENAME_CHR } from '../modules/call_variants/rename_chr'
include { SNPEFF_ANNOTATE_VCF } from '../modules/call_variants/snpeff'

workflow CALLVAR {
    take:
        bam_good

    main:
        chr_name = Channel.value(file(params.chr_name))
        //gff3 = Channel.value(file(params.gff))
        ref = Channel.value(file(params.reference))
        //call_variants = GATK_CALL_VARIANTS(bam_good, ref)//FREEBAYES_CALL_VARIANTS(bam_good, ref)//BCFTOOLS_CALL_VARIANTS(bam_good, ref)
        //call_variants = BCFTOOLS_CALL_VARIANTS(bam_good, ref)
        call_variants = FREEBAYES_CALL_VARIANTS(bam_good, ref)
        vcf = call_variants.other
        bcftools_stats = BCFTOOLS_STATS(vcf)
        other_count = call_variants.other_count
        vcf_renamed = RENAME_CHR(vcf, chr_name)
        vcf_annotated = SNPEFF_ANNOTATE_VCF(vcf_renamed)
        ann = vcf_annotated.ann

    emit:
        vcf
        other_count
        bcftools_stats
        ann
}
