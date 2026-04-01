workflow INPUT_CHECK {
    take:
    samplesheet

    main:
    reads = Channel
        .fromPath(samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def sample = row.sample ?: row.Sample ?: row.id ?: row.ID
            def fastq1 = row.fastq_1 ?: row.R1
            def fastq2 = row.fastq_2 ?: row.R2
            def layout = (row.layout ?: row.Layout ?: '').toString().trim().toUpperCase()
            def isGzFastq = { value ->
                value && value.toString().toLowerCase() ==~ /.*\.(fastq|fq)\.gz$/
            }

            if (!sample) {
                throw new IllegalArgumentException("Samplesheet row is missing sample ID")
            }
            if (!fastq1 && !fastq2) {
                throw new IllegalArgumentException("Samplesheet row for '${sample}' is missing FASTQ paths")
            }
            if (fastq1 && !isGzFastq(fastq1)) {
                throw new IllegalArgumentException("Samplesheet row for '${sample}' has unsupported fastq_1 path '${fastq1}'. Use gzipped FASTQ files (*.fastq.gz or *.fq.gz).")
            }
            if (fastq2 && !isGzFastq(fastq2)) {
                throw new IllegalArgumentException("Samplesheet row for '${sample}' has unsupported fastq_2 path '${fastq2}'. Use gzipped FASTQ files (*.fastq.gz or *.fq.gz).")
            }

            def reads = []
            if (fastq1) {
                reads << file(fastq1, checkIfExists: true)
            }
            if (fastq2) {
                reads << file(fastq2, checkIfExists: true)
            }

            if (layout == 'PAIRED' && reads.size() != 2) {
                throw new IllegalArgumentException("Samplesheet row for '${sample}' is marked PAIRED but does not have exactly two FASTQ files")
            }

            def singleEnd = reads.size() == 1
            tuple([id: sample, single_end: singleEnd], reads)
        }

    emit:
    reads
}
