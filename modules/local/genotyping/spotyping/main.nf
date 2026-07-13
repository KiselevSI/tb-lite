process SPOTYPING {
    tag        "SpoTyping: $sample_name"
    label 'process_low'
    publishDir "${params.outdir}/spotyping/$sample_name", mode: params.mode

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tb-lite/spotyping:2.0' :
        'tb-lite/spotyping:2.0' }"

    input:
        tuple val(sample_name), path(fastq_files)

    output:
        path("$sample_name.*"), emit: other
        path("$sample_name"), emit: code
        path("${sample_name}.tsv"), emit: table
        tuple val("${task.process}"), val('SpoTyping'), eval('SpoTyping.py --version 2>&1 | sed -n "s/^SpoTyping.py[[:space:]]*//p" | head -1'), topic: versions, emit: versions_spotyping

    script:
        def reads = fastq_files instanceof List ? fastq_files : [fastq_files]

        reads = reads.join(" ")

        """
        SpoTyping.py $reads --noQuery -o $sample_name
        awk -F '\\t' -v OFS='\\t' -v S='${sample_name}' \\
            'NR==1 {print S, \$2, \$3}' \\
            ${sample_name} > ${sample_name}.tsv
        """
}
