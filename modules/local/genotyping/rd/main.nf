process RD {
    tag "RD: ${sample_name}"
    label 'process_single'
    publishDir "${params.outdir}/rd/${sample_name}", mode: params.mode

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tb-lite/rd-scanner:2.0' :
        'tb-lite/rd-scanner:2.0' }"

    input:
        tuple val(sample_name), path(bed)
        path db

    output:
        tuple val(sample_name), path("${sample_name}.rd.tsv"), emit: rd
        path("${sample_name}.rd.tsv.errors.tsv")
        tuple val("${task.process}"), val('rd-scanner'), val('2.0'), topic: versions, emit: versions_rd_scanner

    script:
        // Параллелизм даёт Nextflow, поэтому -j 1; --fresh страхует от
        // унаследованного .state.sqlite при повторном запуске задачи.
        """
        rd_scan.py \\
            -i $bed \\
            -k $db \\
            -o ${sample_name}.rd.tsv \\
            -j 1 \\
            --fresh \\
            --no-progress \\
            --log-level INFO

        # rd_scan.py завершается с кодом 0, даже если образец не обработался:
        # причина уходит в errors.tsv, поэтому проверяем его явно.
        if [ "\$(wc -l < ${sample_name}.rd.tsv.errors.tsv)" -gt 1 ]; then
            cat ${sample_name}.rd.tsv.errors.tsv >&2
            exit 1
        fi
        """
}
