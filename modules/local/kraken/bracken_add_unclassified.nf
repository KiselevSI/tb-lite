process BRACKEN_ADD_UNCLASSIFIED {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'ubuntu:22.04' }"

    input:
    tuple val(meta), path(bracken_tsv), path(kraken_report)

    output:
    tuple val(meta), path("*.add_Unclassified.tsv"), emit: tsv
    path "versions.yml"                             , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    add_unclassified.sh -o ${prefix}.bracken.add_Unclassified.tsv ${bracken_tsv} ${kraken_report}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -1 | sed 's/.*version //' | sed 's/ .*//')
    END_VERSIONS
    """
}
