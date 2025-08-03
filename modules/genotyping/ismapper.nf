process ISMAPPER {
    tag        "is6110: $sample_name"
    publishDir "${params.outdir}/is6110/paired", mode: params.mode

    input:
        tuple val(sample_name), path(fastq_files)
        path is6110
        path ref_gbk

    output:
        path("${sample_name}/*")

    script:

        def reads = fastq_files instanceof List ? fastq_files : [fastq_files]
        reads = reads.join(" ")

        """
        ismap --reads $reads --queries $is6110 --reference $ref_gbk --bam --t ${task.cpus} --output_dir $sample_name
        """
}