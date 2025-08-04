process RENAME_CHR {
    tag "run_rename_chromosome_merge_vcf"
    label 'small_mem'
    publishDir "${params.outdir}/TABLE_ANNOTATION", mode: params.mode
    

    input:
        path(vcf)
        path chromosome_name
        

    output:
        path("merged_vcf.renamed_chromosome.vcf")

    script:

        """
        bcftools annotate --rename-chrs $chromosome_name $vcf -o merged_vcf.renamed_chromosome.vcf
        """
}