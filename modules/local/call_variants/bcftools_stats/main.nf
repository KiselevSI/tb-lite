process BCFTOOLS_STATS{
    tag "run_bcftools_stats: ${sample_name}"
    label 'small_mem'
    label 'solo_cpu'
    //publishDir "${params.outdir}/stats/bcftools", mode: params.mode

    input:
        tuple val(sample_name), path(vcf), path(vcf_csi)

    output:
        path("${sample_name}.bcftools.txt")

    script:
        """
        bcftools stats $vcf > ${sample_name}.bcftools.txt
        """
}