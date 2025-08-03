/*
 * filter.nf : TBmix + map/samtools-метрики → список good_samples
 * IN : bam_all  – tuple val(sample_name), path(bam)
 *      reference – path(ref.fa)
 *      levels    – path(levels.tsv)   (загружаем как value-канал)
 * OUT: good_samples – val(sample_name)
 */
include { TB_MIX } from '../modules/filter/tb_mix'
include { RUN_MAP_STATS } from '../modules/filter/map_stats'
include { RUN_SAMTOOLS_STATS } from '../modules/filter/samtools_stats'
workflow FILTER {
    take:
        trimmed
        bam_all
    main:
        ref = Channel.value(file(params.reference))
        levels = Channel.value(file(params.levels))
        tbmix       = TB_MIX(bam_all, ref, levels)

        map_stats   = RUN_MAP_STATS(bam_all, ref)
        samtools_stats    = RUN_SAMTOOLS_STATS(bam_all)

        align_pct = samtools_stats.rmp
        .map { sample_name, rmp_file ->
            tuple( sample_name, rmp_file.text.trim() as Float )
        }
        // ------------------ 2) среднее покрытие ---------------------
        mean_cov = map_stats.mean
            .map { sample_name, cov_file ->
                tuple( sample_name, cov_file.text.trim() as Float )
            }
        // ------------------ 3) медианное покрытие -------------------
        median_cov = map_stats.median
            .map { sample_name, cov_file ->
                tuple( sample_name, cov_file.text.trim() as Float )
            }
        bad_metrics = align_pct
            .join(mean_cov).join(median_cov)
            .filter { id, pct, cov, median ->
                //pct < params.min_align_pct || cov < params.min_mean_cov
                median < params.min_median
            }
            .map { id, pct, cov, median ->
                "${id}\t${pct}\t${cov}\t${median}"
            }

        bad_header = Channel.value("sample_id\treads_mapped_pct\tmean_coverage\tmedian_coverage")

        bad_with_header = bad_header.mix(bad_metrics)

        bad_with_header.collectFile(
            name:     'bad_reads_low_coverage.txt',
            storeDir: params.outdir,
            newLine:  true,
            sort:     false   // сортировать теперь не нужно (мы уже включили header)
        )

        good_samples = median_cov                        // (id, med)
            .filter { id, med -> med >= params.min_median }
            .map    { id, med -> id }                    // val(id)


        trimmed_good  = trimmed.join(good_samples)

        bam_good      = bam_all.join(good_samples)

        paired_reads = trimmed_good.filter { _id, files -> files.size() == 2 }
            .join( good_samples)

        wgs_metrics = map_stats.wgs
        align_metrics = map_stats.align
        samtools_stat = samtools_stats.stats
        samtools_flagstat = samtools_stats.flagstat

    emit:
        trimmed_good
        paired_reads
        bam_good
        tbmix
        wgs_metrics
        align_metrics
        samtools_stat
        samtools_flagstat
}
