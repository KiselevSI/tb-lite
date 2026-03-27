process SNPEFF_ANNOTATE_VCF {
    tag        "SNPEFF_ANNOTATE_VCF: ${sample_name}"
    label 'big_mem'
    label 'solo_cpu'
    publishDir "${params.outdir}/annotate_vcf/${sample_name}", mode: params.mode

    input:
        tuple val(sample_name), path(vcf_renamed)

    output:
        tuple val(sample_name), path("${sample_name}.annotated.vcf.gz"), emit: ann

    script:
        """
        # 1) Аннотация в plain VCF
        snpEff \
          -noLog -noStats -nodownload -no-downstream -no-upstream -no-utr \
          -v ${params.snpeff_db} \
          ${vcf_renamed} \
          > ${sample_name}.annotated.vcf

        # 2) BGZF-сжатие
        bgzip -@ ${task.cpus} -c ${sample_name}.annotated.vcf > ${sample_name}.annotated.vcf.gz
        rm ${sample_name}.annotated.vcf

        # 3) Tabix-индекс
        tabix -p vcf ${sample_name}.annotated.vcf.gz

        """
}
