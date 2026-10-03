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
    path spotyping_logs
    path tblg_table
    path drugs
    path rd
    path gbk
    path spoldb4

    output:
    path "filter.tbmix.tsv"
    path "tbmix.total.tsv"
    path "tblg.total.tsv"
    path "drug_resist.xlsx"
    path "drug_resist_and_uncertain.xlsx"
    path "general.tsv"
    path "spotyping.total.tsv"
    path "spotyping.full.tsv"
    path "spoligo_spacer_counts.tsv"
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
    write_list spotyping_logs.list $spotyping_logs
    write_list tblg_table.list $tblg_table
    write_list drugs.list $drugs
    write_list rd.list $rd

    # Образцы, дошедшие до вызова вариантов: BCFTOOLS_STATS в CALLVAR работает
    # только по bam_good. Метрики покрытия считаются раньше фильтра, поэтому в
    # general.tsv и tbmix.total.tsv иначе попадают отбракованные образцы — в
    # базе они дали бы карточки без линии, сполиготипа, устойчивости и SNP.
    sed 's|.*/||; s|\\.bcftools_stats\\.txt\$||' bcftools_stats.list \\
        | LC_ALL=C sort -u > passed_samples.txt

    python ${projectDir}/bin/build_metrics_table.py \\
        --wgs-list wgs_metrics.list \\
        --bcftools-list bcftools_stats.list \\
        --flagstat-list samtools_flagstat.list \\
        -o general.unfiltered.tsv \\
        --round

    python ${projectDir}/bin/filter_table_by_samples.py \\
        -i general.unfiltered.tsv -s passed_samples.txt --id-column ID -o general.tsv

    python ${projectDir}/bin/concat_tables.py --input-list tbmix.list --keep-header -o tbmix.unfiltered.tsv
    python ${projectDir}/bin/filter_table_by_samples.py \\
        -i tbmix.unfiltered.tsv -s passed_samples.txt --id-column Sample -o tbmix.total.tsv

    python ${projectDir}/bin/concat_tables.py --input-list spotyping.list --prepend-line "Sample\tSpolBin\tSpol8" -o spotyping.total.tsv
    python ${projectDir}/bin/concat_tables.py --input-list tblg_table.list --keep-header -o tblg.total.tsv

    # filter.tbmix.tsv строится из уже отфильтрованного tbmix.total.tsv.
    python ${projectDir}/bin/filter_tbmix.py -f tbmix.total.tsv -t tblg.total.tsv -o filter.tbmix.tsv

    python ${projectDir}/bin/profiler_parser.py --input-list drugs.list -g $gbk --flat-header -o drug_resist.xlsx

    python ${projectDir}/bin/profiler_parser.py --input-list drugs.list -g $gbk --flat-header --include-other-uncertain -o drug_resist_and_uncertain.xlsx

    # Полная spoligo-таблица для сайта: 43 спейсера, min/rmin и SIT/clade (SpolDB4).
    # spotyping.total.tsv остаётся трёхколоночным — его формат ждут импортёры сайта.
    python ${projectDir}/bin/build_spoligo_table.py \\
        --spoligo-list spotyping.list \\
        --log-list spotyping_logs.list \\
        --spoldb4 $spoldb4 \\
        --count tolerant \\
        --spacers-output spoligo_spacer_counts.tsv \\
        --full-output spotyping.full.tsv

    # rd_scan.py уже пишет per-sample TSV в формате сайта — остаётся склеить.
    python ${projectDir}/bin/concat_tables.py --input-list rd.list --keep-header -o rd.tsv
    """
}
