process TB_PROFILER_DR {
    tag "run_tb_profiler_dr: $sample_name"
    publishDir "${params.outdir}/tb-profiler/drug-resist/$sample_name", mode: params.mode

    input:
        tuple val(sample_name), path(vcf), path(vcf_csi)

    output:
        path("*"), emit: other
        path("results/${sample_name}.results.json"), emit: results

    script:

        """
        tb-profiler profile --vcf $vcf -p $sample_name --txt
        """
}