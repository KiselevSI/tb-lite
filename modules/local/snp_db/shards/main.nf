process SNP_DB_SHARDS {
    tag "${shard_tag}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tb-lite/tb-platform-tables:1.1' :
        'tb-lite/tb-platform-tables:1.1' }"

    input:
    tuple val(shard_tag), path(ann_vcfs)

    output:
    path "sites.${shard_tag}.tsv.gz",   emit: sites
    path "alleles.${shard_tag}.tsv.gz", emit: alleles
    path "vcfrows.${shard_tag}.tsv.gz", emit: vcfrows
    path "samples.${shard_tag}.txt",    emit: samples
    path "counts.${shard_tag}.tsv",     emit: counts

    script:
    """
    : > vcfs.list
    for vcf_path in ${ann_vcfs}; do
        printf '%s\\n' "\$(readlink -f "\$vcf_path")" >> vcfs.list
    done

    python ${projectDir}/bin/vcf_to_snp_shards.py \\
        --file-list vcfs.list \\
        --tag ${shard_tag} \\
        -o . \\
        -j ${task.cpus}
    """
}
