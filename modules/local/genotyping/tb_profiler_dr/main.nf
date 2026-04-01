process TB_PROFILER_DR {
    tag "${sample_name}"
    label 'process_medium'

    publishDir "${params.outdir}/tb-profiler/drug-resist", mode: params.mode

    input:
    tuple val(sample_name), path(bam), path(bam_idx)

    output:
    path("${sample_name}/vcf/"), emit: other
    path("${sample_name}/results/${sample_name}.results.json"), emit: json

    script:
    """
    samtools view -H ${bam} \
        | sed "s/SN:NC_000962\\.3/SN:Chromosome/" \
        | samtools reheader - ${bam} > ${sample_name}.dedup.chromosome.bam

    samtools index ${sample_name}.dedup.chromosome.bam

    tb-profiler profile \
        --threads ${task.cpus} \
        --bam ${sample_name}.dedup.chromosome.bam \
        --prefix ${sample_name} \
        --txt \
        --call_whole_genome

    mkdir -p ${sample_name}
    mv bam ${sample_name}/
    mv results ${sample_name}/
    mv vcf ${sample_name}/
    mv ${sample_name}.dedup.chromosome.bam ${sample_name}/
    mv ${sample_name}.dedup.chromosome.bam.bai ${sample_name}/
    """
}
