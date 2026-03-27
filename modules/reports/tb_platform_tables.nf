process TB_PLATFORM_TABLES {
    tag "Final Table"
    label 'small_mem'
    publishDir "${params.outdir}/tb-platform", mode: params.mode

    input:
    path multiqc
    path tbmix
    path spotyping
    path tblg_table
    path drugs
    path rd
    path gbk

    output:
    path "filter.tbmix.tsv"
    path "dr.xlsx"
    path "dr_other_variants.xlsx"
    path "general.tsv"
    path "spotyping.total.tsv"
    path "rd.tsv"
    

    script:
    """

    awk -F '\\t' '(NR==1) || (FNR>1)' $tbmix  > tbmix.total.tsv
    

    echo -e "Sample\\tSpolBin\\tSpol8" > spotyping.total.tsv

    cat $spotyping  >> spotyping.total.tsv

    awk -F '\\t' '(NR==1) || (FNR>1)' $tblg_table  > tblg.total.tsv

    filter_tbmix.py -f tbmix.total.tsv -t tblg.total.tsv -o filter.tbmix.tsv

    profiler_parser.py -i $drugs -g $gbk --flat-header --also-other-variants -o dr.xlsx

    deletions_to_csv.py -i $rd -o rd.tsv
  
    build_general_table.py --round -m $multiqc -o general.tsv
    
    """
}