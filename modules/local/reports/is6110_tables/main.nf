process IS6110_TABLES {
    tag "IS6110 Tables"
    label 'process_single'
    publishDir "${params.outdir}/Reports/tb-platform/is6110", mode: params.mode

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tb-lite/tb-platform-tables:1.1' :
        'tb-lite/tb-platform-tables:1.1' }"

    input:
    path is6110_tables

    output:
    path "*.tsv"

    script:
    """
    write_list() {
        local output_file="\$1"
        shift
        : > "\$output_file"
        for input_path in "\$@"; do
            [[ -n "\$input_path" ]] && printf '%s\\n' "\$input_path" >> "\$output_file"
        done
        if [[ -s "\$output_file" ]]; then
            LC_ALL=C sort -u "\$output_file" -o "\$output_file"
        fi
    }

    write_list is6110.list $is6110_tables

    python ${projectDir}/bin/build_is6110_tables.py \\
        --input-list is6110.list \\
        --outdir .
    """
}
