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
    def krakenInputs = kraken_combined ? kraken_combined.collect { it.getName() }.join(' ') : ''
    def tableInputs = ['general.tsv', 'spotyping.total.tsv', 'filter.tbmix.tsv', 'tblg.total.tsv', 'tbmix.total.tsv']
    if (kraken_combined) {
        tableInputs << 'kraken.top_hits.tsv'
    }
    def buildFinalArgs = tableInputs.join(' ')
    """
    python ${projectDir}/bin/build_metrics_table.py --wgs $wgs_metrics --bcftools $bcftools_stats --flagstat $samtools_flagstat -o general.tsv --round

    {
        echo "ID"
        for f in $fastp_reports; do
            basename "\$f" .fastp.json
        done | sort -u
    } > all_samples.tsv

    awk -F '\\t' '(NR==1) || (FNR>1)' $tbmix  > tbmix.total.tsv


    echo -e "Sample\\tSpolBin\\tSpol8" > spotyping.total.tsv

    cat $spotyping  >> spotyping.total.tsv

    awk -F '\\t' '(NR==1) || (FNR>1)' $tblg_table  > tblg.total.tsv

    python ${projectDir}/bin/filter_tbmix.py -f tbmix.total.tsv -t tblg.total.tsv -o filter.tbmix.tsv

    python ${projectDir}/bin/profiler_parser.py -i $drugs -g $gbk --include-other-uncertain -o drug_resist.xlsx

    ${ kraken_combined ? "python ${projectDir}/bin/kraken_top_hits.py -i ${krakenInputs} -o kraken.top_hits.tsv" : "" }

    python ${projectDir}/bin/build_final_table.py --base-table all_samples.tsv --dr drug_resist.xlsx -t ${buildFinalArgs} --join left --str-cols Spol8,SpolBin -o FINAL_TABLE.xlsx


    """
}
