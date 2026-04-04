process TB_MIX {
    tag "run_lineage: ${sample_name}"
    label 'process_single'
    publishDir "${params.outdir}/tb-mix", mode: params.mode

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tb-lite/tb-mix:1.0' :
        'tb-lite/tb-mix:1.0' }"

    input:
        tuple val(sample_name), path(bam), path(bai)
        path ref
        path levels

    output:
        path("${sample_name}.mix.tsv")

    script:
        """
        tb_mix.py -i $bam -r $ref -l $levels --mq 30 --bq 20 -f 0.04 -p $sample_name -o ${sample_name}.mix.tsv --mix-low 0.05 --mix-high 0.95

        """
}
