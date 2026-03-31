process SNPEFF_ANN {
    tag        "SNPEFF_ANNOTATE_VCF"
    label 'big_mem'
    label 'solo_cpu'
    publishDir "${params.outdir}/snpeff", mode: params.mode

    input:
        path(vcf)
        path chromosome_name

    output:
        tuple path("cohort.ann.vcf.gz"), path("cohort.ann.vcf.gz.tbi")

    script:
        """
        # 1) Аннотация в plain VCF
        # bcftools annotate --rename-chrs $chromosome_name $vcf -O z -o renamed_chromosome.vcf.gz
        snpEff \
          -noLog -noStats -nodownload -no-downstream -no-upstream -no-utr \
          -v ${params.snpeff_db} \
          $vcf | bcftools view -Oz -o cohort.ann.vcf.gz

        bcftools index -t -f cohort.ann.vcf.gz

        """
}
