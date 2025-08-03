process MULTIQC {
    tag "multiqc"
    publishDir "${params.outdir}/multiqc", mode: params.mode

    input:
    path reports
    path cfg

    output:
    path "*"
    path "tb_multiqc_report_data/multiqc_data.json", emit: report

    script:
    """
    multiqc --config $cfg --title 'TB-Lite QC' --filename tb_multiqc_report .
    """
}