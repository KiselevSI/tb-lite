
process POST_PROCESS_TABLE {
    tag        "SNPEFF_table"
    label 'big_mem'
    label 'solo_cpu'
    publishDir "${params.outdir}/snpeff", mode: params.mode

    input:
        path tsv
        path feature_table

    output:
        path "FINAL_ANNOTATION_TABLE.tsv", emit: clean

    script: 
    """
    awk -v OFS="\t" 'BEGIN{FS="\t"} {for(i=1;i<=NF;i++) if(\$i==".") \$i=""; print}' "${tsv}" > FINAL_ANNOTATION_TABLE_1.tsv
    add_name_strand.py -i FINAL_ANNOTATION_TABLE_1.tsv -t $feature_table -o FINAL_ANNOTATION_TABLE_2.tsv
    dedup_ann_columns.py -i FINAL_ANNOTATION_TABLE_2.tsv -o FINAL_ANNOTATION_TABLE.tsv
    """

}
