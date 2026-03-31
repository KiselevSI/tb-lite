process TBLG {
    tag        "tblg: $sample_name"
    label 'small_mem'
    label 'solo_cpu'
    publishDir "${params.outdir}/lineage", mode: params.mode

    input:
        tuple val(sample_name), path(vcf), path(vcf_csi)

    output:
        path("${sample_name}.lg.tsv")

    script:
        
        """
        tblg $vcf -o ${sample_name}.lg.tsv
        """
}