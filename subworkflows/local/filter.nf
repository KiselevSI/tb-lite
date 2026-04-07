/*
 * filter.nf : TBmix + map/samtools-метрики → список good_samples
 * IN : bam_all  – tuple val(sample_name), path(bam)
 *      reference – path(ref.fa)
 *      levels    – path(levels.tsv)   (загружаем как value-канал)
 * OUT: good_samples – val(sample_name)
 */
include { TB_MIX } from '../../modules/local/filter/tb_mix/main'
include { PICARD_COLLECTWGSMETRICS } from '../../modules/nf-core/picard/collectwgsmetrics/main'
include { PICARD_COLLECTALIGNMENTSUMMARYMETRICS } from '../../modules/nf-core/picard/collectalignmentsummarymetrics/main'
include { SAMTOOLS_STATS } from '../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/samtools/flagstat/main'
workflow FILTER {
    take:
        trimmed
        bam_all
        ref_fai
    main:
        ref = Channel.value(file(params.reference))
        levels = Channel.value(file(params.levels))
        tbmix       = TB_MIX(bam_all, ref, levels)

        ref_meta = [id: file(params.reference).baseName]
        ref_tuple = Channel.value([ref_meta, file(params.reference)])

        ch_bam_meta = bam_all.map { sample_name, bam, bai ->
            [[id: sample_name], bam, bai]
        }
        ch_bam_meta_no_index = bam_all.map { sample_name, bam, bai ->
            [[id: sample_name], bam]
        }

        PICARD_COLLECTWGSMETRICS(
            ch_bam_meta,
            ref_tuple,
            ref_fai,
            []
        )

        PICARD_COLLECTALIGNMENTSUMMARYMETRICS(
            ch_bam_meta_no_index,
            ref_tuple
        )

        ch_fasta_fai = ref_fai.map { meta, fai ->
            [meta, file(params.reference), fai]
        }

        SAMTOOLS_STATS(ch_bam_meta, ch_fasta_fai)
        SAMTOOLS_FLAGSTAT(ch_bam_meta)

        wgs_metrics = PICARD_COLLECTWGSMETRICS.out.metrics
            .map { meta, metrics -> metrics }
        wgs_metrics_with_id = PICARD_COLLECTWGSMETRICS.out.metrics
            .map { meta, metrics -> tuple(meta.id, metrics) }
        align_metrics = PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.metrics
            .map { meta, metrics -> metrics }
        samtools_stat = SAMTOOLS_STATS.out.stats
            .map { meta, stats -> stats }
        samtools_flagstat = SAMTOOLS_FLAGSTAT.out.flagstat
            .map { meta, flagstat -> flagstat }
        samtools_flagstat_with_id = SAMTOOLS_FLAGSTAT.out.flagstat
            .map { meta, flagstat -> tuple(meta.id, flagstat) }

        align_pct = samtools_flagstat_with_id
            .map { sample_name, flagstat_file ->
                def text = flagstat_file.text
                def match = text =~ /mapped \(([0-9]+(?:\.[0-9]+)?)%/
                def pct = match ? match[0][1] as Float : 0.0f
                tuple(sample_name, pct)
            }
        // ------------------ 2) среднее покрытие ---------------------
        mean_cov = wgs_metrics_with_id
            .map { sample_name, metrics_file ->
                def data_line = metrics_file.readLines().findAll { it && !it.startsWith('#') }[1]
                def parts = data_line.split('\t')
                tuple(sample_name, parts[1] as Float)
            }
        // ------------------ 3) медианное покрытие -------------------
        median_cov = wgs_metrics_with_id
            .map { sample_name, metrics_file ->
                def data_line = metrics_file.readLines().findAll { it && !it.startsWith('#') }[1]
                def parts = data_line.split('\t')
                tuple(sample_name, parts[3] as Float)
            }
        combined_metrics = align_pct
            .join(mean_cov)
            .join(median_cov)

        bad_metrics = combined_metrics
            .filter { id, pct, cov, median ->
                pct <= params.min_align_pct || median < params.min_median
            }
            .map { id, pct, cov, median ->
                "${id}\t${pct}\t${cov}\t${median}"
            }

        bad_header = Channel.value("sample_id\treads_mapped_pct\tmean_coverage\tmedian_coverage")

        bad_with_header = bad_header.mix(bad_metrics)
        bad_reads_name = params.batch_tag
            ? "bad_reads_low_coverage.${params.batch_tag}.txt"
            : 'bad_reads_low_coverage.txt'
        bad_reads_dir = params.batch_tag
            ? "${params.outdir}/batch_reports/filter"
            : "${params.outdir}/Reports/general"

        bad_with_header.collectFile(
            name:     bad_reads_name,
            storeDir: bad_reads_dir,
            newLine:  true,
            sort:     false   // сортировать теперь не нужно (мы уже включили header)
        )

        good_samples = combined_metrics
            .filter { id, pct, cov, median ->
                pct >= params.min_align_pct && median >= params.min_median
            }
            .map { id, pct, cov, median -> id }


        trimmed_good  = trimmed.join(good_samples)

        bam_good      = bam_all.join(good_samples)

    emit:
        trimmed_good
        bam_good
        tbmix
        wgs_metrics
        align_metrics
        samtools_stat
        samtools_flagstat
}
