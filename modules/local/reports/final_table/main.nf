process FINAL_TABLE {
    tag "Final Table"
    label 'process_single'
    publishDir "${params.outdir}/Reports/general",
        mode: params.mode,
        saveAs: { filename -> ['FINAL_TABLE.xlsx', 'drug_resist.xlsx'].contains(filename) ? filename : null }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tb-lite/build-table:1.1' :
        'tb-lite/build-table:1.1' }"

    input:
    path wgs_metrics
    path fastp_reports
    path bcftools_stats
    path samtools_flagstat
    path tbmix
    path spotyping
    path tblg_table
    path drugs
    path gbk
    path kraken_combined

    output:
    path "*"


    script:
    def tableInputs = ['general.tsv', 'spotyping.total.tsv', 'filter.tbmix.tsv', 'tblg.total.tsv', 'tbmix.total.tsv']
    if (kraken_combined) {
        tableInputs << 'kraken.top_hits.tsv'
    }
    def buildFinalArgs = tableInputs.join(' ')
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
    write_list fastp_reports.list $fastp_reports
    write_list bcftools_stats.list $bcftools_stats
    write_list samtools_flagstat.list $samtools_flagstat
    write_list tbmix.list $tbmix
    write_list spotyping.list $spotyping
    write_list tblg_table.list $tblg_table
    write_list drugs.list $drugs
    write_list kraken_combined.list $kraken_combined

    python ${projectDir}/bin/build_metrics_table.py \\
        --wgs-list wgs_metrics.list \\
        --bcftools-list bcftools_stats.list \\
        --flagstat-list samtools_flagstat.list \\
        -o general.tsv \\
        --round

    {
        echo "ID"
        while IFS= read -r input_path; do
            base_name="\${input_path##*/}"
            printf '%s\\n' "\${base_name%.fastp.json}"
        done < fastp_reports.list | LC_ALL=C sort -u
    } > all_samples.tsv

    python ${projectDir}/bin/concat_tables.py --input-list tbmix.list --keep-header -o tbmix.total.tsv
    python ${projectDir}/bin/concat_tables.py --input-list spotyping.list --prepend-line "Sample\tSpolBin\tSpol8" -o spotyping.total.tsv
    python ${projectDir}/bin/concat_tables.py --input-list tblg_table.list --keep-header -o tblg.total.tsv

    python ${projectDir}/bin/filter_tbmix.py -f tbmix.total.tsv -t tblg.total.tsv -o filter.tbmix.tsv

    python ${projectDir}/bin/profiler_parser.py --input-list drugs.list -g $gbk --include-other-uncertain -o drug_resist.xlsx

    ${ kraken_combined ? "python ${projectDir}/bin/kraken_top_hits.py --input-list kraken_combined.list -o kraken.top_hits.tsv" : "" }

    python ${projectDir}/bin/build_final_table.py --base-table all_samples.tsv --dr drug_resist.xlsx -t ${buildFinalArgs} --join left --str-cols Spol8,SpolBin -o FINAL_TABLE.xlsx
    """
}
