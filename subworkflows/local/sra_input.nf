include { CUSTOM_SRATOOLSNCBISETTINGS } from '../../modules/nf-core/custom/sratoolsncbisettings/main'
include { SRATOOLS_PREFETCH } from '../../modules/nf-core/sratools/prefetch/main'
include { SRATOOLS_FASTERQDUMP } from '../../modules/nf-core/sratools/fasterqdump/main'

process SRA_DETECT_LAYOUT {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::sra-tools=3.2.1 conda-forge::pigz=2.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-tools:3.2.1--h4304569_1' :
        'quay.io/biocontainers/sra-tools:3.2.1--h4304569_1' }"

    input:
    tuple val(meta), path(sra)

    output:
    tuple val(meta), path(sra), path('layout.txt'), optional: true, emit: supported
    tuple val(meta), path('unsupported_layout.tsv'), optional: true, emit: unsupported

    script:
    """
    if [[ -f "${sra}" && "${sra}" == *.sra ]]; then
        sra_file="${sra}"
    else
        sra_file=\$(find -L ${sra} -type f -name '*.sra' | head -n 1)
    fi

    if [[ -z "\$sra_file" ]]; then
        echo "ERROR: prefetched SRA file for ${meta.id} not found under ${sra}" >&2
        exit 1
    fi

    mkdir -p preview

    fasterq-dump \\
        --threads 1 \\
        --split-files \\
        --skip-technical \\
        --outdir preview \\
        "\$sra_file"

    shopt -s nullglob
    fastqs=(preview/*.fastq)
    count=\${#fastqs[@]}

    if [[ "\$count" -eq 1 ]]; then
        echo single > layout.txt
    elif [[ "\$count" -eq 2 ]]; then
        echo paired > layout.txt
    else
        printf 'sample_id\\treason\\tfastq_count\\n%s\\tunsupported_layout\\t%s\\n' "${meta.id}" "\$count" > unsupported_layout.tsv
        echo "WARN: ${meta.id} produced \$count FASTQ files during layout detection; skipping accession." >&2
    fi
    """
}

workflow SRA_INPUT {
    take:
    ids_file

    main:
    accessions = Channel
        .fromPath(ids_file)
        .splitText()
        .map { it.trim() }
        .filter { it && !it.startsWith('#') }
        .unique()
        .map { accession -> [[id: accession], accession] }

    accession_ids = Channel
        .fromPath(ids_file)
        .splitText()
        .map { it.trim() }
        .filter { it && !it.startsWith('#') }
        .unique()

    CUSTOM_SRATOOLSNCBISETTINGS(accession_ids.collect())

    ncbi_settings = CUSTOM_SRATOOLSNCBISETTINGS.out.ncbi_settings
    certificate = Channel.value([])

    SRATOOLS_PREFETCH(accessions, ncbi_settings, certificate)
    SRA_DETECT_LAYOUT(SRATOOLS_PREFETCH.out.sra)

    unsupported_name = params.batch_tag
        ? "unsupported_sra_layout.${params.batch_tag}.txt"
        : 'unsupported_sra_layout.txt'
    unsupported_dir = params.batch_tag
        ? "${params.outdir}/batch_reports/filter"
        : "${params.outdir}/Reports/general"

    unsupported_header = Channel.value("sample_id\treason\tfastq_count")
    unsupported_rows = SRA_DETECT_LAYOUT.out.unsupported
        .flatMap { meta, report ->
            report.text.readLines().drop(1).findAll { it }
        }

    unsupported_header
        .mix(unsupported_rows)
        .collectFile(
            name: unsupported_name,
            storeDir: unsupported_dir,
            newLine: true,
            sort: false
        )

    fasterq_input = SRA_DETECT_LAYOUT.out.supported.map { meta, sra, layout ->
        def singleEnd = layout.text.trim() == 'single'
        [meta + [single_end: singleEnd], sra]
    }

    SRATOOLS_FASTERQDUMP(fasterq_input, ncbi_settings, certificate)

    reads = SRATOOLS_FASTERQDUMP.out.reads.map { meta, fastqs ->
        def files = fastqs instanceof List ? fastqs.sort { it.name } : [fastqs]
        tuple([id: meta.id, single_end: meta.single_end], files)
    }

    emit:
    reads
}
