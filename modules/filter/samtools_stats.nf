process RUN_SAMTOOLS_STATS{
    tag "run_samtools_stats: ${sample_name}"
    label 'small_mem'
    label 'solo_cpu'
    publishDir "${params.outdir}/stats/samtools/${sample_name}", mode: params.mode

    input:
        tuple val(sample_name), path(bam), path(bai)

    output:
        path("${sample_name}.samtools.txt"), emit: stats
        path("${sample_name}.flagstat.txt"),  emit: flagstat
        tuple val(sample_name), path("${sample_name}.rmp.txt"), emit: rmp

    script:
        """
        samtools stats ${bam} > ${sample_name}.samtools.txt
        samtools flagstat ${bam} > ${sample_name}.flagstat.txt

        perc=\$(grep -oP '\\(\\K[0-9]+\\.[0-9]+(?=%)' \
            ${sample_name}.flagstat.txt | head -n1)
        if [ -z \"\$perc\" ]; then perc="0.0"; fi
        echo \"\$perc\" > ${sample_name}.rmp.txt

        """
}