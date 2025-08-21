process SNPEFF_ANNOTATE_VCF {
    tag        "SNPEFF_ANNOTATE_VCF: ${sample_name}"
    label 'big_mem'
    label 'solo_cpu'
    publishDir "${params.outdir}/annotate_vcf/${sample_name}", mode: params.mode

    input:
        tuple val(sample_name), path(vcf_renamed)        

    output:
        tuple val(sample_name), path("${sample_name}.annotated.vcf.gz"), emit: ann
        tuple val(sample_name), path("${sample_name}.ann.tsv"), emit: ann_table

    script:

        """
        snpEff ann -noLog -noStats -no-downstream -no-upstream -no-utr -v Mycobacterium_tuberculosis_h37rv ${vcf_renamed} | bgzip -c > ${sample_name}.annotated.vcf.gz
        
        bcftools index -f ${sample_name}.annotated.vcf.gz
        
        bcftools query -f '%POS\\t%REF\\t%ALT\\t%QUAL\\t%INFO/ANN\\n' \
            ${sample_name}.annotated.vcf.gz > ${sample_name}.ann.tsv

        """
}