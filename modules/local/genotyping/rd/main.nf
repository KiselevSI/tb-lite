process RD {
    tag "RD: ${sample_name}"
    label 'process_single'
    publishDir "${params.outdir}/rd/${sample_name}", mode: params.mode

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tb-lite/rd-scanner:1.0' :
        'tb-lite/rd-scanner:1.0' }"

    input:
        tuple val(sample_name), path(bed)
        path db

    output:
        tuple val(sample_name), path("${sample_name}.novel_rd.tsv"), emit: rd
        path("${sample_name}.known_rd.tsv")
        tuple val("${task.process}"), val('rd-scanner'), val('1.0'), topic: versions, emit: versions_rd_scanner

    script:
        """
        rd.py $bed \
        -k $db \
        -n ${sample_name}.novel_rd.tsv \
        -o ${sample_name}.known_rd.tsv
        """
}
