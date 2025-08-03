process ANNOTATION_TABLE {
    tag "run_rename_chromosome_merge_vcf"
    publishDir "${params.outdir}/TABLE_ANNOTATION", mode: params.mode
    

    input:
        path(vcfs)
        path vcf_ann
        path feature_table
        

    output:
        path("ANNOTATION_TABLE.xlsx")

    script:

        """
        vcf2table.py -v $vcfs -a $vcf_ann -t $feature_table -o ANNOTATION_TABLE.xlsx
        """
}