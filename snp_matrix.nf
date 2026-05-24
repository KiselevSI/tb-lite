#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SAMTOOLS_FAIDX } from './modules/nf-core/samtools/faidx/main'
include { ANN_TABLE } from './subworkflows/local/ann_table'

params.vcf_input = null

process PREPARE_SNP_MATRIX_VCF {
    tag "${sample_name}"
    label 'process_low'

    conda "${projectDir}/modules/nf-core/bcftools/view/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
        tuple val(sample_name), path(input_vcf)

    output:
        tuple val(sample_name), path("${sample_name}.snp_matrix.vcf.gz"), path("${sample_name}.snp_matrix.vcf.gz.tbi"), path("${sample_name}.var_count.txt")

    script:
    """
    set -euo pipefail

    sample_count=\$(bcftools query -l "${input_vcf}" | wc -l | tr -d '[:space:]')
    if [ "\${sample_count}" -ne 1 ]; then
        echo "Input VCF for '${sample_name}' must contain exactly one sample column, found \${sample_count}: ${input_vcf}" >&2
        exit 1
    fi

    printf "%s\\n" "${sample_name}" > sample_name.txt
    bcftools view -Oz -o normalized.vcf.gz "${input_vcf}"
    bcftools reheader -s sample_name.txt -o "${sample_name}.snp_matrix.vcf.gz" normalized.vcf.gz
    bcftools index -t "${sample_name}.snp_matrix.vcf.gz"
    bcftools view -H "${sample_name}.snp_matrix.vcf.gz" | wc -l | tr -d '[:space:]' > "${sample_name}.var_count.txt"
    """
}

if (params.help) {
    log.info """
    ============================
    TB-Lite SNP Matrix Workflow
    ============================

    Usage:
      nextflow run snp_matrix.nf --vcf_input <vcf_samples.csv> [options]

    Required:
      --vcf_input           CSV with columns: sample,vcf

    Input CSV example:
      sample,vcf
      sample1,/data/sample1.vcf
      sample2,/data/sample2.vcf.gz

    Options:
      --outdir              Output directory [default: ./results2]
      --reference           Reference FASTA [default: assets/h37rv.fa]
      --snpeff_db           SnpEff database name [default: h37rv_custom]
      --snpeff_data_dir     Path to SnpEff data directory [default: assets/SNPEFF_ANNOTATION/data]
      --h37rv_feature_table Feature table for final SNP matrix post-processing
      --mode                publishDir mode: copy or link [default: copy]
      --help                Show this help message

    Notes:
      Input VCFs can be .vcf or .vcf.gz and must contain exactly one sample column.
      Per-sample annotated VCF is not required; annotation is applied after cohort merge.
    """.stripIndent()
    System.exit(0)
}

if (!params.vcf_input) {
    log.error "ERROR: set --vcf_input to a CSV file with columns: sample,vcf. Use --help for usage information."
    System.exit(1)
}

workflow SNP_MATRIX_FROM_VCF {
    main:
        vcf_records = Channel
            .fromPath(params.vcf_input, checkIfExists: true)
            .splitCsv(header: true, sep: ',')
            .map { row ->
                def sample = (row.sample ?: row.Sample ?: row.id ?: row.ID)?.toString()?.trim()
                def vcf = (row.vcf ?: row.VCF)?.toString()?.trim()

                if (!sample) {
                    throw new IllegalArgumentException("VCF input row is missing sample ID")
                }
                if (!(sample ==~ /[A-Za-z0-9._-]+/)) {
                    throw new IllegalArgumentException("Sample ID '${sample}' contains unsupported characters. Use only letters, numbers, '.', '_' and '-'.")
                }
                if (!vcf) {
                    throw new IllegalArgumentException("VCF input row for '${sample}' is missing VCF path")
                }

                tuple(sample, file(vcf, checkIfExists: true))
            }

        PREPARE_SNP_MATRIX_VCF(vcf_records)

        ref_meta = [id: file(params.reference).baseName]
        ref_ch = Channel.value([ref_meta, file(params.reference, checkIfExists: true)])
        ch_faidx_input = ref_ch.map { meta, fasta -> [meta, fasta, []] }
        SAMTOOLS_FAIDX(ch_faidx_input, false)

        ANN_TABLE(PREPARE_SNP_MATRIX_VCF.out, SAMTOOLS_FAIDX.out.fai)

    emit:
        clean = ANN_TABLE.out.clean
}

workflow {
    SNP_MATRIX_FROM_VCF()
}
