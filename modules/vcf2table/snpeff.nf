process ANNOTATION_VCFS {
    tag "run_annotation_merge_vcf"
    label 'small_mem'
    publishDir "${params.outdir}/TABLE_ANNOTATION", mode: params.mode

    input:
        path merged_vcf
        

    output:
        path("merged.annotation.vcf")

    script:
        """
        snpEff ann -noLog -noStats -no-downstream -no-upstream -no-utr -v Mycobacterium_tuberculosis_h37rv $merged_vcf > merged.annotation.vcf
        """
}