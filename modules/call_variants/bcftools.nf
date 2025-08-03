process BCFTOOLS_CALL_VARIANTS {
    tag "call_variants: ${sample_name}"
    publishDir "${params.outdir}/vcf/${sample_name}", mode: params.mode

    input:
        tuple val(sample_name), path(bam), path(bam_idx)
        path ref
        

    output:
        tuple val(sample_name), path("${sample_name}.vcf.gz"), path("${sample_name}.vcf.gz.csi"), emit:other
        path("${sample_name}.vcf.gz"), emit:vcfs

    script:
        """
        bcftools mpileup --threads ${task.cpus} \
            --min-MQ 30 --ignore-overlaps --max-depth 3000 \
            -f ${ref} ${bam} -Ou | \
        bcftools call --threads ${task.cpus} --multiallelic-caller \
            --ploidy 1 --variants-only -Ou - | \
        bcftools view --threads ${task.cpus} \
            --include 'QUAL>20 && DP>10' -Oz -o ${sample_name}.vcf.gz
        bcftools index ${sample_name}.vcf.gz
        """
}