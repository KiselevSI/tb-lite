process FINAL_TABLE {
    tag "Final Table"
    publishDir "${params.outdir}", mode: params.mode

    input:
    path multiqc
    path tbmix
    path spotyping
    path tblg_table
    path drugs

    output:
    path "*"
    

    script:
    """
    awk -F '\\t' '(NR==1) || (FNR>1)' $tbmix  > tbmix.total.tsv

    echo -e "Sample\\tSpolBin\\tSpol8" > spotyping.total.tsv

    cat $spotyping  >> spotyping.total.tsv

    awk -F '\\t' '(NR==1) || (FNR>1)' $tblg_table  > tblg.total.tsv

    filter_tbmix.py -f tbmix.total.tsv -t tblg.total.tsv -o filter.tbmix.tsv

    dr_parser.py -i $drugs -o dr.xlsx
  
    build_final_table.py -m $multiqc --dr dr.xlsx -t spotyping.total.tsv filter.tbmix.tsv tblg.total.tsv tbmix.total.tsv --str-cols Spol8,SpolBin -o FINAL_TABLE.xlsx

    
    """
}