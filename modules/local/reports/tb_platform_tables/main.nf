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
    path "drug_resist.xlsx"
    path "drug_resist_and_uncertain.xlsx"
    path "general.tsv"
    path "spotyping.total.tsv"
    path "rd.tsv"


    script:
    """
    python ${projectDir}/bin/build_metrics_table.py --wgs $wgs_metrics --bcftools $bcftools_stats --flagstat $samtools_flagstat -o general.tsv --round

    awk -F '\\t' '(NR==1) || (FNR>1)' $tbmix  > tbmix.total.tsv


    echo -e "Sample\\tSpolBin\\tSpol8" > spotyping.total.tsv

    cat $spotyping  >> spotyping.total.tsv

    awk -F '\\t' '(NR==1) || (FNR>1)' $tblg_table  > tblg.total.tsv

    python ${projectDir}/bin/filter_tbmix.py -f tbmix.total.tsv -t tblg.total.tsv -o filter.tbmix.tsv

    python ${projectDir}/bin/profiler_parser.py -i $drugs -g $gbk --flat-header -o drug_resist.xlsx

    python ${projectDir}/bin/profiler_parser.py -i $drugs -g $gbk --flat-header --include-other-uncertain -o drug_resist_and_uncertain.xlsx

    python ${projectDir}/bin/deletions_to_csv.py -i $rd -o rd.tsv

    """
}
