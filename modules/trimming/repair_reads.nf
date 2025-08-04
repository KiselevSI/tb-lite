process REPAIR_READS {
    tag "repair: $sample_name"
    publishDir "${params.outdir}/repaired/$sample_name", mode: params.mode

    input:
    tuple val(sample_name), path(fastq_files)

    output:
    tuple val(sample_name), path("${sample_name}_*.repaired.fastq.gz"), emit: repaired_reads

    script:
    // Предполагаем только парные риды
    def (r1, r2) = fastq_files
    """
    repair.sh \
        in1=$r1 in2=$r2 \
        out1=${sample_name}_1.repaired.fastq.gz \
        out2=${sample_name}_2.repaired.fastq.gz \
        outs=${sample_name}.orphans.fastq.gz \
        repair tossbrokenreads
    """
}