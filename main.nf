
/* === Подключаем все sub-workflow-ы ============================ */
include { TRIMMING }      from './subworkflows/trimming'
include { QC       }      from './subworkflows/qc'
include { MAPPING  }      from './subworkflows/mapping'
include { FILTER   }      from './subworkflows/filter'
include { CALLVAR  }      from './subworkflows/call_varints'
include { GENOTYPE }      from './subworkflows/genotyping'
include { REPORTS  }      from './subworkflows/reports'
include { ANN_TABLE  }      from './subworkflows/ann_table'

/* === Top-level workflow ====================================== */
workflow {

    sample_count = Math.max(
        file(params.samples).toFile().readLines().findAll { it.trim() }.size() - 1,
        0
    )

    reads = Channel
      .fromPath(params.samples)           // файл с таблицей
      .splitCsv(header:true, sep:',')    // Map на каждую строку
      .map { row ->                      // row.Sample, row.R1, row.R2, row.Layout
          def sample  = row.Sample
          def layout  = row.Layout?.toUpperCase()
          def r1      = file(row.R1)      // преобразуем в Path
          def fq_list = (layout == 'PAIRED' && row.R2)
                        ? [ r1, file(row.R2) ]   // список из двух файлов
                        : [ r1 ]                 // список из одного
          tuple(sample, fq_list)          // итоговый элемент канала
      }

    /* 2. TRIMMING --------------------------------------------------------- */
    trim = TRIMMING( reads )                               // → TRIMMING.trimmed_reads

    /* 3. FASTQC ----------------------------------------------------------- */
    qc = QC( trim.trimmed_reads )                        // → QC.fastqc_reports

    /* 4. ALIGNMENT -------------------------------------------------------- */
    maps = MAPPING( trim.trimmed_reads)                         // → MAPPING.bam

    /* 5. FILTERING + МЕТРИКИ --------------------------------------------- */
    filt = FILTER( trim.trimmed_reads,
            maps.bam) 


    /* 6. VARIANT CALLING -------------------------------------------------- */
    callvar = CALLVAR( filt.bam_good)

    /* 7. GENOTYPING ------------------------------------------------------- */


    gen = GENOTYPE(filt.trimmed_good,
              callvar.vcf,
              filt.bam_good )
    /* 8. ОТЧЁТЫ ----------------------------------------------------------- */
    if (!params.skip_reports) {
        REPORTS(
            filt.wgs_metrics,         // wgs_metrics
            filt.align_metrics,       // align_metrics
            qc.reports,          // fastqc_reports
            callvar.bcftools_stats,     // bcftools_stats
            filt.samtools_stat,       // samtools_stats
            filt.samtools_flagstat,   // samtools_flagstat
            filt.tbmix,               // tbmix
            gen.spol_table,        // spotyping
            gen.tblg_table,        // tblg_table  (должен быть emit’ом GENOTYPE)
            gen.dr,
            gen.del
        )
    }

    if (sample_count > 1 & !params.skip_reports) {
        ANN_TABLE(
            callvar.other_count
        )
    } else {
        log.info "Skipping ANN_TABLE: input sample count is ${sample_count}, need more than 1"
    }


}
