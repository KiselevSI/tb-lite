
process POST_PROCESS_TABLE {
    tag        "SNPEFF_table"
    label 'process_low'
    publishDir "${params.outdir}/Reports/snp_matrix", mode: params.mode

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tb-lite/ann-table:1.0' :
        'tb-lite/ann-table:1.0' }"

    input:
        path tsv
        path feature_table

    output:
        path "FINAL_ANNOTATION_TABLE.tsv", emit: clean

    script:
    """
    awk -v OFS="\t" 'BEGIN{FS="\t"} {for(i=1;i<=NF;i++) if(\$i==".") \$i=""; print}' "${tsv}" > FINAL_ANNOTATION_TABLE_1.tsv
    python ${projectDir}/bin/add_name_strand.py -i FINAL_ANNOTATION_TABLE_1.tsv -t $feature_table -o FINAL_ANNOTATION_TABLE_2.tsv
    python ${projectDir}/bin/dedup_ann_columns.py -i FINAL_ANNOTATION_TABLE_2.tsv -o FINAL_ANNOTATION_TABLE.tsv
    """

}
