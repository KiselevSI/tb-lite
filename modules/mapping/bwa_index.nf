process BWA_INDEX {
    tag "index reference: $ref"
    label 'small_mem'
    publishDir("${params.outdir}/ref/", mode: params.mode)
    
    input:
    path ref
    
    output:
    path "${ref}.*", emit: index
    path ref, emit: ref_file
    script:

    """
        bwa index $ref
    """


}