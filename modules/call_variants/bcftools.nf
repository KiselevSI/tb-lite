process BCFTOOLS_CALL_VARIANTS {
    tag "call_variants: ${sample_name}"
    label 'big_mem'
    label 'multi_cpu'
    publishDir "${params.outdir}/vcf/${sample_name}", mode: params.mode

    input:
        tuple val(sample_name), path(bam), path(bam_idx)
        path ref
        

    output:
        tuple val(sample_name), path("${sample_name}.vcf.gz"), path("${sample_name}.vcf.gz.csi"), emit:other
        path("${sample_name}.vcf.gz"), emit:vcfs

    script:
        """
        bcftools mpileup \
            --threads ${task.cpus} \
            -f ${ref} \
            -q 20 -Q 20 \
            --max-depth 10000 \
            -a AD,DP \
            ${bam} -Ou | \
        bcftools call \
            --threads ${task.cpus} \
            -m --ploidy 1 -v -Ou | \
        bcftools norm \
            -f ${ref} -m -both -c w -Ou | \
        bcftools view \
            --threads ${task.cpus} \
            -i 'QUAL>=20 && FMT/DP>=10' \
            -Oz -o ${sample_name}.vcf.gz

        bcftools index -f ${sample_name}.vcf.gz
        """
}