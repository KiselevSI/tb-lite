process FREEBAYES_CALL_VARIANTS {
    tag "call_variants: ${sample_name}"
    label 'big_mem'
    label 'multi_cpu'
    publishDir "${params.outdir}/vcf/${sample_name}", mode: params.mode

    input:
        tuple val(sample_name), path(bam), path(bam_idx)
        path ref
        

    output:
        tuple val(sample_name), path("${sample_name}.vcf.gz"), path("${sample_name}.vcf.gz.tbi"), emit:other
        tuple val(sample_name), path("${sample_name}.vcf.gz"), path("${sample_name}.vcf.gz.tbi"), path("${sample_name}.var_count.txt"), emit:other_count

    script:
        """
        freebayes \
            --fasta-reference $ref \
            --ploidy 1 \
            --min-mapping-quality 20 \
            --min-base-quality 20 \
            --min-coverage 10 \
            --min-alternate-fraction 0.80 --max-complex-gap 3 \
            -b $bam | bcftools norm -f $ref -Oz --threads ${task.cpus} | \
            bcftools filter -e 'QUAL<10 || INFO/DP<10' -Oz --threads ${task.cpus} -o ${sample_name}.vcf.gz

            bcftools index -t -f ${sample_name}.vcf.gz

            zgrep -cv '#' ${sample_name}.vcf.gz > ${sample_name}.var_count.txt
        """
}