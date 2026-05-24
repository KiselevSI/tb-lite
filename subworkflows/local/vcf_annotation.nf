include { BCFTOOLS_ANNOTATE as VCF_ANNOTATION_BCFTOOLS_ANNOTATE } from '../../modules/nf-core/bcftools/annotate/main'
include { SNPEFF_SNPEFF as VCF_ANNOTATION_SNPEFF } from '../../modules/nf-core/snpeff/snpeff/main'

process PREPARE_VCF_ANNOTATION_INPUT {
    tag "${meta.id}"
    label 'process_low'

    conda "${projectDir}/modules/nf-core/bcftools/view/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(input_vcf)

    output:
    tuple val(meta), path("${meta.id}.vcf_annotation_input.vcf.gz"), path("${meta.id}.vcf_annotation_input.vcf.gz.tbi"), emit: vcf

    script:
    """
    set -euo pipefail

    sample_count=\$(bcftools query -l "${input_vcf}" | wc -l | tr -d '[:space:]')
    if [[ "\${sample_count}" -ne 1 ]]; then
        echo "Input VCF for '${meta.id}' must contain exactly one sample column, found \${sample_count}: ${input_vcf}" >&2
        exit 1
    fi

    printf "%s\\n" "${meta.id}" > sample_name.txt
    bcftools view -Oz -o normalized.vcf.gz "${input_vcf}"
    bcftools reheader -s sample_name.txt -o "${meta.id}.vcf_annotation_input.vcf.gz" normalized.vcf.gz
    bcftools index -t "${meta.id}.vcf_annotation_input.vcf.gz"
    """

    stub:
    """
    echo "" | gzip > ${meta.id}.vcf_annotation_input.vcf.gz
    touch ${meta.id}.vcf_annotation_input.vcf.gz.tbi
    """
}

workflow VCF_ANNOTATION {
    take:
    vcf_list

    main:
    vcf_records = Channel
        .fromPath(vcf_list, checkIfExists: true)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def sample = (row.sample ?: row.Sample ?: row.id ?: row.ID)?.toString()?.trim()
            def vcf = (row.vcf ?: row.VCF)?.toString()?.trim()

            if (!sample) {
                throw new IllegalArgumentException("VCF list row is missing sample ID")
            }
            if (!(sample ==~ /[A-Za-z0-9._-]+/)) {
                throw new IllegalArgumentException("Sample ID '${sample}' contains unsupported characters. Use only letters, numbers, '.', '_' and '-'.")
            }
            if (!vcf) {
                throw new IllegalArgumentException("VCF list row for '${sample}' is missing VCF path")
            }
            if (!(vcf.toLowerCase() ==~ /.*\.vcf(\.gz)?$/)) {
                throw new IllegalArgumentException("VCF list row for '${sample}' has unsupported VCF path '${vcf}'. Use *.vcf or *.vcf.gz.")
            }

            tuple([id: sample], file(vcf, checkIfExists: true))
        }

    PREPARE_VCF_ANNOTATION_INPUT(vcf_records)

    snpeff_cache = Channel.value([[id: 'snpeff_cache'], file(params.snpeff_data_dir, checkIfExists: true)])

    if (params.snpeff_db == "Mycobacterium_tuberculosis_h37rv") {
        chr_name = file(params.chr_name, checkIfExists: true)
        renamed_input = PREPARE_VCF_ANNOTATION_INPUT.out.vcf.map { meta, vcf, tbi ->
            [meta, vcf, tbi, [], [], [], [], chr_name]
        }
        VCF_ANNOTATION_BCFTOOLS_ANNOTATE(renamed_input)
        snpeff_input = VCF_ANNOTATION_BCFTOOLS_ANNOTATE.out.vcf
    } else {
        snpeff_input = PREPARE_VCF_ANNOTATION_INPUT.out.vcf.map { meta, vcf, tbi ->
            [meta, vcf]
        }
    }

    VCF_ANNOTATION_SNPEFF(snpeff_input, params.snpeff_db, snpeff_cache)

    ann = VCF_ANNOTATION_SNPEFF.out.vcf.map { meta, ann_vcf ->
        tuple(meta.id, ann_vcf)
    }

    emit:
    ann
}
