process TB_PROFILER_DR {
    tag "run_tb_profiler_dr: $sample_name"
    label 'big_mem'
    label 'multi_cpu'
    publishDir "${params.outdir}/tb-profiler/drug-resist/$sample_name", mode: params.mode

    input:
        tuple val(sample_name), path(bam), path(bam_idx) 

    output:
        path("*"), emit: other
        path("results/${sample_name}.results.json"), emit: results

    script:

        """
        tb-profiler profile -t ${task.cpus} --bam $bam -p $sample_name --txt --call_whole_genome
        """
}