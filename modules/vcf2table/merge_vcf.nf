process MERGE_VCF {
    tag "run_merge_vcf"
    publishDir "${params.outdir}/TABLE_ANNOTATION", mode: params.mode

    input:
        path vcfs
        

    output:
        path("merged.vcf")

    script:
        """
        merge_vcf.py -v $vcfs -o merged.vcf
        """
}