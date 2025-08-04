process RUN_MAP_STATS {
    tag "map_stats: ${sample_name}"
    label 'medium_mem'
    label 'solo_cpu'
    publishDir "${params.outdir}/stats/picard/${sample_name}", mode: params.mode

    input:
        tuple val(sample_name), path(bam), path(bai)
        path ref

    output:
        path("${sample_name}.wgs_metrics.txt"), emit: wgs
        path("${sample_name}.alignment_metrics.txt"), emit: align
        tuple val(sample_name), path("${sample_name}.mean_coverage.txt"), emit: mean
        tuple val(sample_name), path("${sample_name}.median_coverage.txt"), emit: median

    script:
        """
        java -jar /usr/picard/picard.jar \
        CollectWgsMetrics COVERAGE_CAP=100000 \
        I=$bam \
        O=${sample_name}.wgs_metrics.txt \
        R=$ref COUNT_UNPAIRED=true

        java -jar /usr/picard/picard.jar \
          CollectAlignmentSummaryMetrics \
          R=$ref \
          I=$bam \
          O=${sample_name}.alignment_metrics.txt
        grep -v "^#" ${sample_name}.wgs_metrics.txt | awk 'NR==3 {print \$2}'> ${sample_name}.mean_coverage.txt
        grep -v "^#" ${sample_name}.wgs_metrics.txt | awk 'NR==3 {print \$4}' > ${sample_name}.median_coverage.txt
        """

}