include { CUSTOM_SRATOOLSNCBISETTINGS } from '../../modules/nf-core/custom/sratoolsncbisettings/main'
include { SRATOOLS_PREFETCH } from '../../modules/nf-core/sratools/prefetch/main'
include { SRATOOLS_FASTERQDUMP } from '../../modules/nf-core/sratools/fasterqdump/main'

process SRA_DETECT_LAYOUT {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(sra)

    output:
    tuple val(meta), path("${meta.id}"), path('layout.txt')

    script:
    """
    sra_file=\$(find ${sra} -type f -name '*.sra' | head -n 1)

    if [[ -z "\$sra_file" ]]; then
        echo "ERROR: prefetched SRA file for ${meta.id} not found under ${sra}" >&2
        exit 1
    fi

    mkdir -p preview

    fasterq-dump \\
        --threads 1 \\
        --split-files \\
        --skip-technical \\
        -X 1 \\
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
        echo "ERROR: ${meta.id} produced \$count FASTQ files during layout detection; expected 1 or 2." >&2
        exit 1
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

    fasterq_input = SRA_DETECT_LAYOUT.out.map { meta, sra, layout ->
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
