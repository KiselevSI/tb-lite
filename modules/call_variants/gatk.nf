process GATK_CALL_VARIANTS {
    tag "gatk_call_variants: ${sample_name}"
    label 'big_mem'
    label 'multi_cpu'
    publishDir "${params.outdir}/vcf_gatk/${sample_name}", mode: params.mode

    input:
        tuple val(sample_name), path(bam), path(bam_idx)
        path ref

    output:
        tuple val(sample_name), path("${sample_name}.vcf.gz"), path("${sample_name}.vcf.gz.tbi"), emit: other
        tuple path("${sample_name}.vcf.gz"), path("${sample_name}.vcf.gz.tbi"), emit: vcfs

    script:
    """
    samtools faidx "${ref}"
    gatk CreateSequenceDictionary -R ${ref} -O ${ref.baseName}.dict


    # Вызов вариантов (haploid)
    gatk HaplotypeCaller \
        -R "${ref}" \
        -I "${bam}" \
        -O "${sample_name}.raw.vcf.gz" \
        --sample-ploidy 1 \
        --min-base-quality-score 20 \
        --standard-min-confidence-threshold-for-calling 30 \
        --native-pair-hmm-threads ${task.cpus}

    # Жёсткие фильтры GATK по аннотациям качества
    gatk VariantFiltration \
        -R "${ref}" \
        -V "${sample_name}.raw.vcf.gz" \
        -O "${sample_name}.flt.vcf.gz" \
        --filter-name "LowQD"      --filter-expression "QD < 2.0" \
        --filter-name "HighFS"     --filter-expression "FS > 60.0" \
        --filter-name "HighSOR"    --filter-expression "SOR > 3.0" \
        --filter-name "LowMQ"      --filter-expression "MQ < 40.0" \
        --filter-name "LowMQRank"  --filter-expression "MQRankSum < -12.5" \
        --filter-name "LowRPRS"    --filter-expression "ReadPosRankSum < -8.0"

    # Оставить только PASS
    gatk SelectVariants -V "${sample_name}.flt.vcf.gz" --exclude-filtered -O "${sample_name}.pass.vcf.gz"

    # Нормализация + близкие к вашим пороги (QUAL, DP) и доля ALT ≈ 0.80 по FORMAT/AD
    bcftools norm -f "${ref}" -m -both -Ou "${sample_name}.pass.vcf.gz" | \
      bcftools filter -e 'QUAL<10 || INFO/DP<10' \
        -Oz --threads ${task.cpus} -o "${sample_name}.vcf.gz"

    bcftools index -t -f "${sample_name}.vcf.gz"
    """
}
