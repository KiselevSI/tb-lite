process RD {
    tag "RD: ${sample_name}"
    publishDir "${params.outdir}/rd/${sample_name}", mode: params.mode

    input:
        tuple val(sample_name), path(bed)
        path db

    output:
        tuple val(sample_name), path("${sample_name}.novel_rd.tsv"),
        path("${sample_name}.known_rd.tsv")

    script:
        """
        rd.py $bed \
        -k $db \
        -n ${sample_name}.novel_rd.tsv \
        -o ${sample_name}.known_rd.tsv
        """
}