include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { SRA_INPUT }   from '../subworkflows/local/sra_input'
include { TRIMMING }    from '../subworkflows/local/trimming'
include { QC }          from '../subworkflows/local/qc'
include { MAPPING }     from '../subworkflows/local/mapping'
include { FILTER }      from '../subworkflows/local/filter'
include { CALLVAR }     from '../subworkflows/local/call_variants'
include { GENOTYPE }    from '../subworkflows/local/genotyping'
include { KRAKEN }      from '../subworkflows/local/kraken'
include { ANN_TABLE }   from '../subworkflows/local/ann_table'
include { REPORTS }     from '../subworkflows/local/reports'
include { VCF_ANNOTATION } from '../subworkflows/local/vcf_annotation'
include { VERSIONS }     from '../subworkflows/local/versions'

workflow TBLITE {
    main:
    resolved_input = params.input ?: params.samples
    skip_multiqc = params.skip_multiqc || params.skip_reports
    skip_final_reports = params.skip_final_reports || params.skip_reports

    ch_reads = Channel.empty()
    ch_fastqc = Channel.empty()
    ch_kraken_multiqc = Channel.empty()
    ch_kraken_combined = Channel.empty()
    ch_extra_versions = Channel.empty()

    if (params.vcf_annotation_only) {
        VCF_ANNOTATION(params.vcf_list)
    } else {
        if (resolved_input) {
            INPUT_CHECK(resolved_input)
            ch_reads = INPUT_CHECK.out.reads
        } else {
            SRA_INPUT(params.sra_ids)
            ch_reads = SRA_INPUT.out.reads
        }

        trim = TRIMMING(ch_reads)

        if (!params.skip_qc) {
            qc = QC(trim.trimmed_reads)
            ch_fastqc = qc.reports
        }

        maps = MAPPING(trim.trimmed_reads_legacy)
        filt = FILTER(trim.trimmed_reads_legacy, maps.bam, maps.ref_fai)

        if (!params.skip_kraken && params.kraken2_db) {
            kraken_reads = trim.trimmed_reads.map { meta, reads ->
                [meta + [single_end: reads.size() == 1], reads]
            }
            kraken = KRAKEN(kraken_reads)
            ch_kraken_multiqc = kraken.multiqc
            ch_kraken_combined = kraken.combined
        }

        callvar = CALLVAR(filt.bam_good, maps.ref_fai)
        gen = GENOTYPE(filt.trimmed_good, callvar.vcf, filt.bam_good)
        if (!params.skip_snp_matrix) {
            ANN_TABLE(callvar.other_count, maps.ref_fai)
        }

        if (!skip_multiqc || !skip_final_reports) {
            reports = REPORTS(
                filt.wgs_metrics,
                filt.align_metrics,
                trim.fastp_json,
                ch_fastqc,
                callvar.bcftools_stats,
                filt.samtools_stat,
                filt.samtools_flagstat,
                ch_kraken_multiqc,
                ch_kraken_combined,
                filt.tbmix,
                gen.spol_table,
                gen.spol_log,
                gen.tblg_table,
                gen.dr,
                gen.del,
                gen.is6110
            )
            ch_extra_versions = reports.versions
        }
    }

    VERSIONS(ch_extra_versions)
}
