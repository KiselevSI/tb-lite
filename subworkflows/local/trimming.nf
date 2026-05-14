/*
 * trimming.nf : nf-core/fastp
 * IN : tuple val(meta), path(reads)
 * OUT: trimmed_reads – tuple val(meta), path(*.fastp.fastq.gz)
 */
include { FASTP } from '../../modules/nf-core/fastp/main'

process VALIDATE_RAW_FASTQ {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::gzip=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/55/556474e164daf5a5e218cd5d497681dcba0645047cf24698f88e3e078eacbd09/data' :
        'community.wave.seqera.io/library/fastp:1.1.0--08aa7c5662a30d57' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('valid/*'), optional: true, emit: valid
    tuple val(meta), path('invalid_fastq.tsv'), optional: true, emit: invalid

    script:
    def read_files = reads instanceof List ? reads : [reads]
    def read_args = read_files.collect { "\"${it}\"" }.join(' ')
    """
    read_files=(${read_args})
    bad_entries=()

    for fq in "\${read_files[@]}"; do
        if [[ ! -s "\$fq" ]]; then
            bad_entries+=("\$(basename "\$fq"):empty_file")
        elif ! gzip -t "\$fq" >/dev/null 2>&1; then
            bad_entries+=("\$(basename "\$fq"):invalid_gzip")
        elif [[ "\$(gzip -cd "\$fq" | head -c 1 | wc -c)" -eq 0 ]]; then
            bad_entries+=("\$(basename "\$fq"):empty_fastq")
        fi
    done

    if [[ "\${#bad_entries[@]}" -gt 0 ]]; then
        files=\$(IFS=,; echo "\${bad_entries[*]}")
        printf 'sample_id\\treason\\tfiles\\n%s\\traw_invalid_fastq\\t%s\\n' "${meta.id}" "\$files" > invalid_fastq.tsv
        echo "WARN: ${meta.id} has invalid raw FASTQ input: \$files; skipping sample." >&2
        exit 0
    fi

    mkdir -p valid
    for fq in "\${read_files[@]}"; do
        ln -s "../\$fq" "valid/\$(basename "\$fq")"
    done
    """
}

process VALIDATE_TRIMMED_FASTQ {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::gzip=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/55/556474e164daf5a5e218cd5d497681dcba0645047cf24698f88e3e078eacbd09/data' :
        'community.wave.seqera.io/library/fastp:1.1.0--08aa7c5662a30d57' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('valid/*'), optional: true, emit: valid
    tuple val(meta), path('invalid_fastq.tsv'), optional: true, emit: invalid

    script:
    def read_files = reads instanceof List ? reads : [reads]
    def read_args = read_files.collect { "\"${it}\"" }.join(' ')
    """
    read_files=(${read_args})
    bad_entries=()

    for fq in "\${read_files[@]}"; do
        if [[ ! -s "\$fq" ]]; then
            bad_entries+=("\$(basename "\$fq"):empty_file")
        elif ! gzip -t "\$fq" >/dev/null 2>&1; then
            bad_entries+=("\$(basename "\$fq"):invalid_gzip")
        elif [[ "\$(gzip -cd "\$fq" | head -c 1 | wc -c)" -eq 0 ]]; then
            bad_entries+=("\$(basename "\$fq"):empty_fastq")
        fi
    done

    if [[ "\${#bad_entries[@]}" -gt 0 ]]; then
        files=\$(IFS=,; echo "\${bad_entries[*]}")
        printf 'sample_id\\treason\\tfiles\\n%s\\ttrimmed_invalid_fastq\\t%s\\n' "${meta.id}" "\$files" > invalid_fastq.tsv
        echo "WARN: ${meta.id} has invalid trimmed FASTQ output: \$files; skipping sample." >&2
        exit 0
    fi

    mkdir -p valid
    for fq in "\${read_files[@]}"; do
        ln -s "../\$fq" "valid/\$(basename "\$fq")"
    done
    """
}

workflow TRIMMING {
    take:
        raw_reads

    main:
        VALIDATE_RAW_FASTQ(raw_reads)

        ch_fastp_input = VALIDATE_RAW_FASTQ.out.valid.map { meta, reads -> [meta, reads, []] }

        FASTP(ch_fastp_input, false, false, false)

        VALIDATE_TRIMMED_FASTQ(FASTP.out.reads)

        trimmed_reads = VALIDATE_TRIMMED_FASTQ.out.valid
        trimmed_reads_legacy = trimmed_reads.map { meta, reads -> tuple(meta.id, reads) }
        fastp_json = FASTP.out.json.map { meta, json -> json }

        invalid_name = params.batch_tag
            ? "bad_reads_invalid_fastq.${params.batch_tag}.txt"
            : 'bad_reads_invalid_fastq.txt'
        invalid_dir = params.batch_tag
            ? "${params.outdir}/batch_reports/filter"
            : "${params.outdir}/Reports/general"

        invalid_header = Channel.value("sample_id\treason\tfiles")
        invalid_rows = VALIDATE_RAW_FASTQ.out.invalid
            .mix(VALIDATE_TRIMMED_FASTQ.out.invalid)
            .flatMap { meta, report ->
                report.text.readLines().drop(1).findAll { it }
            }

        invalid_header
            .mix(invalid_rows)
            .collectFile(
                name: invalid_name,
                storeDir: invalid_dir,
                newLine: true,
                sort: false
            )

    emit:
        trimmed_reads
        trimmed_reads_legacy
        fastp_json
}
