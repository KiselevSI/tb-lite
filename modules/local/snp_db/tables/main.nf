process SNP_DB_TABLES {
    tag "SNP DB tables"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tb-lite/tb-platform-tables:1.1' :
        'tb-lite/tb-platform-tables:1.1' }"

    input:
    path shards

    output:
    path "snp_sites.tsv.gz",          emit: sites
    path "sample_snp_alleles.tsv.gz", emit: alleles
    path "vcf_table.tsv.gz",          emit: vcf_table
    path "samples.txt",               emit: samples
    path "snp_db_manifest.tsv",       emit: manifest
    path "import_snp.sql",            emit: import_sql

    script:
    // Шарды уже разложены Nextflow'ом в рабочем каталоге задачи, поэтому
    // --shard-dir . Итоговые имена (snp_sites.*, samples.txt) под маски
    // шардов (sites.*.tsv.gz, samples.*.txt) не попадают.
    def publish_dir = file(params.outdir).toAbsolutePath().toString() + '/Reports/tb-platform/snp'
    """
    bash ${projectDir}/bin/merge_snp_shards.sh \\
        --shard-dir . \\
        --out-dir . \\
        --threads ${task.cpus}

    # Путь в SQL — это место публикации, откуда файлы будет читать psql.
    python ${projectDir}/bin/write_import_snp_sql.py \\
        --data-dir "${publish_dir}" \\
        -o import_snp.sql
    """
}
