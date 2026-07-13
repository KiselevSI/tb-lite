workflow VERSIONS {
    take:
    extra_versions

    main:
    pipeline_versions = Channel.of(
        tuple('TB-Lite', workflow.manifest.version),
        tuple('Nextflow', workflow.nextflow.version)
    )

    topic_versions = Channel.topic('versions')
        .map { _process_name, program, version ->
            tuple(program, version)
        }

    topic_versions
        .mix(extra_versions)
        .mix(pipeline_versions)
        .map { program, version ->
            def normalizedProgram = program?.toString()?.trim()?.replaceAll(/[\t\r\n]+/, ' ')
            def normalizedVersion = version?.toString()?.trim()?.replaceAll(/[\t\r\n]+/, ' ')
            tuple(normalizedProgram, normalizedVersion)
        }
        .filter { program, version -> program && version }
        .unique()
        .collect(flat: false)
        .map { versions ->
            def rows = versions
                .collect { program, version -> "${program}\t${version}" }
                .unique()
                .sort()
            (['program\tversion'] + rows).join('\n')
        }
        .collectFile(
            name: 'versions.txt',
            storeDir: params.outdir,
            newLine: true
        )
}
