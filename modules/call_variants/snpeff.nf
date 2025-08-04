process SNPEFF_ANNOTATE_VCF {
    tag        "drug_resist: $sample_name"
    label 'small_mem'
    label 'solo_cpu'
    publishDir "${params.outdir}/annotate_vcf", mode: params.mode

    input:
        tuple val(sample_name), path(vcf_renamed)        

    output:
        tuple val(sample_name), path("${sample_name}.annotated.vcf.gz")

    script:

        """
        snpEff ann -noLog -noStats -no-downstream -no-upstream -no-utr -v Mycobacterium_tuberculosis_h37rv $vcf_renamed | bgzip -c > ${sample_name}.annotated.vcf.gz
        """
}