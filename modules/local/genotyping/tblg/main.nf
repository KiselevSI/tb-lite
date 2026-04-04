process TBLG {
    tag        "tblg: $sample_name"
    label 'process_single'
    publishDir "${params.outdir}/lineage", mode: params.mode

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tb-lite/tblg:1.0' :
        'tb-lite/tblg:1.0' }"

    input:
        tuple val(sample_name), path(vcf), path(vcf_csi)

    output:
        path("${sample_name}.lg.tsv")

    script:

        """
        tblg $vcf -o ${sample_name}.lg.tsv
        """
}
