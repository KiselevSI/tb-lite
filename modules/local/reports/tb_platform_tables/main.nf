process TB_PLATFORM_TABLES {
    tag "Final Table"
    label 'process_single'
    publishDir "${params.outdir}/Reports/tb-platform", mode: params.mode

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tb-lite/tb-platform-tables:1.1' :
        'tb-lite/tb-platform-tables:1.1' }"

    input:
    path wgs_metrics
    path bcftools_stats
    path samtools_flagstat
    path tbmix
    path spotyping
    path tblg_table
    path drugs
    path rd
    path gbk

    output:
    path "filter.tbmix.tsv"
    path "tbmix.total.tsv"
    path "drug_resist.xlsx"
    path "drug_resist_and_uncertain.xlsx"
    path "general.tsv"
    path "spotyping.total.tsv"
    path "rd.tsv"


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

    write_list wgs_metrics.list $wgs_metrics
    write_list bcftools_stats.list $bcftools_stats
    write_list samtools_flagstat.list $samtools_flagstat
    write_list tbmix.list $tbmix
    write_list spotyping.list $spotyping
    write_list tblg_table.list $tblg_table
    write_list drugs.list $drugs
    write_list rd.list $rd

    python ${projectDir}/bin/build_metrics_table.py \\
        --wgs-list wgs_metrics.list \\
        --bcftools-list bcftools_stats.list \\
        --flagstat-list samtools_flagstat.list \\
        -o general.tsv \\
        --round

    python ${projectDir}/bin/concat_tables.py --input-list tbmix.list --keep-header -o tbmix.total.tsv
    python ${projectDir}/bin/concat_tables.py --input-list spotyping.list --prepend-line "Sample\tSpolBin\tSpol8" -o spotyping.total.tsv
    python ${projectDir}/bin/concat_tables.py --input-list tblg_table.list --keep-header -o tblg.total.tsv

    python ${projectDir}/bin/filter_tbmix.py -f tbmix.total.tsv -t tblg.total.tsv -o filter.tbmix.tsv

    python ${projectDir}/bin/profiler_parser.py --input-list drugs.list -g $gbk --flat-header -o drug_resist.xlsx

    python ${projectDir}/bin/profiler_parser.py --input-list drugs.list -g $gbk --flat-header --include-other-uncertain -o drug_resist_and_uncertain.xlsx

    python ${projectDir}/bin/deletions_to_csv.py --input-list rd.list -o rd.tsv
    """
}
