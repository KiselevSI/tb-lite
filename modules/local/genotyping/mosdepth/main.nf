process MOSDEPTH {
    tag "mosdepth: ${sample_name}"
    label 'small_mem'
    label 'multi_cpu'
    //publishDir "${params.outdir}/stats/mosdepth/${sample_name}", mode: params.mode

    input:
        tuple val(sample_name), path(bam), path(bam_idx)        

    output:
        tuple path("${sample_name}.mosdepth.global.dist.txt"), path("${sample_name}.per-base.bed.gz.csi")

        tuple val(sample_name), path("${sample_name}.mosdepth.summary.txt"), emit: cov

        tuple val(sample_name), path("${sample_name}.per-base.bed.gz"),  emit: bed

        tuple val(sample_name), path("${sample_name}.median.txt"), emit: median



    script:
        """
        mosdepth -t ${task.cpus} $sample_name $bam

        sort -k2,2nr ${sample_name}.mosdepth.global.dist.txt | awk '\$3>=0.5 {print \$2; exit}' > ${sample_name}.median.txt

        """
}