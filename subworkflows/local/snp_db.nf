include { SNP_DB_SHARDS } from '../../modules/local/snp_db/shards/main'
include { SNP_DB_SHARDS as REBUILD_SHARDS } from '../../modules/local/snp_db/shards/main'
include { SNP_DB_TABLES } from '../../modules/local/snp_db/tables/main'

// Детерминированный тег по составу пачки: перезапуск перезаписывает тот же
// файл, а не создаёт второй с дублями. Имена образцов в пачках не пересекаются,
// поэтому первого имени и размера достаточно для уникальности.
def snpDbChunkTag(List paths, String prefix) {
    def names = paths.collect { it.name }.sort()
    // tokenize('.')[0] вместо регулярки: в slashy-строке Groovy символ '$'
    // перед '/' даёт ошибку компиляции.
    return "${prefix}_${names.first().tokenize('.')[0]}_${names.size()}"
}

workflow SNP_DB {
    take:
        ann_vcfs   // канал путей к *.annotated.ann.vcf

    main:
        // run_batches.sh уже передаёт --batch_tag batch_${i}; без него тег
        // выводится из состава батча.
        tagged = ann_vcfs
            .collect()
            .map { paths ->
                tuple(params.batch_tag ? params.batch_tag as String
                                       : snpDbChunkTag(paths, 'batch'), paths)
            }

        SNP_DB_SHARDS(tagged)

    emit:
        sites   = SNP_DB_SHARDS.out.sites
        alleles = SNP_DB_SHARDS.out.alleles
        vcfrows = SNP_DB_SHARDS.out.vcfrows
}

workflow SNP_DB_FINAL {
    main:
        // Наличие шардов — это состояние файловой системы на момент запуска,
        // а не поток данных. Решаем в Groovy: так каждый канал потребляется
        // ровно один раз и ветки не конфликтуют.
        def shard_dir = file("${params.outdir}/snp_db/shards")
        def existing_shards = shard_dir.exists() ? file("${shard_dir}/sites.*.tsv.gz") : []

        if (existing_shards) {
            log.info("SNP DB: найдено ${existing_shards.size()} шардов в ${shard_dir}")
            shards = Channel
                .fromPath("${shard_dir}/{sites,alleles,vcfrows,samples,counts}.*")
                .collect()
        }
        else {
            log.info("SNP DB: шардов нет, собираю их из ${params.outdir}/annotate_vcf")
            chunks = Channel
                .fromPath("${params.outdir}/annotate_vcf/*/*.annotated.ann.vcf", checkIfExists: true)
                .collate(500)
                .map { chunk -> tuple(snpDbChunkTag(chunk, 'rebuilt'), chunk) }

            REBUILD_SHARDS(chunks)

            shards = REBUILD_SHARDS.out.sites
                .mix(
                    REBUILD_SHARDS.out.alleles,
                    REBUILD_SHARDS.out.vcfrows,
                    REBUILD_SHARDS.out.samples,
                    REBUILD_SHARDS.out.counts
                )
                .collect()
        }

        SNP_DB_TABLES(shards)

    emit:
        manifest = SNP_DB_TABLES.out.manifest
}
