class WorkflowMain {

    public static void initialise(workflow, params, log) {
        def effectiveInput = params.input ?: params.samples
        def effectiveMultiqcConfig = params.multiqc_config ?: params.multiqc
        def skipMultiqc = params.skip_multiqc || params.skip_reports
        def skipFinalReports = params.skip_final_reports || params.skip_reports
        def skipSnpMatrix = params.skip_snp_matrix
        def effectiveProfile = (!workflow.profile || workflow.profile == 'standard') ? 'docker (default runtime)' : workflow.profile

        if (params.samples && !params.input) {
            log.warn "Parameter --samples is deprecated; use --input instead."
        }

        if (params.multiqc && !params.multiqc_config) {
            log.warn "Parameter --multiqc is deprecated; use --multiqc_config instead."
        }

        if (params.skip_reports) {
            log.warn "Parameter --skip_reports is deprecated; use --skip_multiqc and --skip_final_reports instead."
        }

        // Print help message
        if (params.help) {
            log.info """
            ======================
            TB-Lite Pipeline v${workflow.manifest.version}
            ======================

            Usage:
              nextflow run main.nf --input <samplesheet.csv> [options]
              nextflow run main.nf --sra_ids <ids.txt> [options]
              nextflow run main.nf --vcf_annotation_only --vcf_list <vcf_samples.csv> [options]

            Required:
              WGS mode: exactly one of:
                --input       CSV file with sample information (sample,fastq_1,fastq_2; FASTQ must be *.fastq.gz or *.fq.gz)
                --sra_ids     Text file with one SRA accession per line
              VCF annotation-only mode:
                --vcf_annotation_only
                --vcf_list    CSV file with columns sample,vcf

            Options:
              --outdir              Output directory [default: ./results]
              --reference           Reference genome FASTA [default: assets/h37rv.fa]
              --snpeff_data_dir     Path to SnpEff data directory [default: assets/SNPEFF_ANNOTATION/data]
              --kraken2_db          Path to first Kraken2 database
              --kraken2_db_label    Optional label for first Kraken2 database
              --kraken2_db_2        Optional path to second Kraken2 database
              --kraken2_db_label_2  Optional label for second Kraken2 database
              --skip_kraken         Skip Kraken2/Bracken branch [default: false]
              --skip_multiqc        Skip MultiQC generation [default: false]
              --skip_final_reports  Skip final TB reports [default: false]
              --skip_snp_matrix     Skip cohort SNP matrix generation [default: false]
              --skip_qc             Skip FastQC [default: false]
              --mode                publishDir mode: copy or link [default: copy]
              --help                Show this help message

            Profiles:
              default              Docker runtime is enabled even without -profile
              -profile docker      Docker runtime
              -profile singularity Singularity / Apptainer runtime
              -profile conda       Conda runtime

            Kubernetes:
              use -profile docker -c k8s.config
            """.stripIndent()
            System.exit(0)
        }

        if (params.vcf_annotation_only) {
            if (!params.vcf_list) {
                log.error "ERROR: set --vcf_list when using --vcf_annotation_only. Use --help for usage information."
                System.exit(1)
            }
            if (effectiveInput || params.sra_ids) {
                log.error "ERROR: do not set --input, --samples, or --sra_ids with --vcf_annotation_only."
                System.exit(1)
            }
        } else {
            if ((effectiveInput && params.sra_ids) || (!effectiveInput && !params.sra_ids)) {
                log.error "ERROR: set exactly one of --input or --sra_ids. Use --help for usage information."
                System.exit(1)
            }
            if (params.vcf_list) {
                log.warn "Parameter --vcf_list is ignored unless --vcf_annotation_only is set."
            }
        }

        // Log pipeline info
        log.info """
        ======================
        TB-Lite Pipeline v${workflow.manifest.version}
        ======================
        input              : ${effectiveInput ?: '(none)'}
        sra_ids            : ${params.sra_ids ?: '(none)'}
        vcf_annotation_only: ${params.vcf_annotation_only}
        vcf_list           : ${params.vcf_list ?: '(none)'}
        outdir             : ${params.outdir}
        reference          : ${params.reference}
        snpeff_data_dir    : ${params.snpeff_data_dir}
        kraken2_db         : ${params.kraken2_db ?: '(disabled)'}
        kraken2_db_label   : ${params.kraken2_db_label ?: '(auto)'}
        kraken2_db_2       : ${params.kraken2_db_2 ?: '(disabled)'}
        kraken2_db_label_2 : ${params.kraken2_db_label_2 ?: '(auto)'}
        skip_qc            : ${params.skip_qc}
        skip_kraken        : ${params.skip_kraken}
        skip_multiqc       : ${skipMultiqc}
        skip_final_reports : ${skipFinalReports}
        skip_snp_matrix    : ${skipSnpMatrix}
        multiqc_config     : ${effectiveMultiqcConfig ?: '(default)'}
        profile            : ${effectiveProfile}
        ======================
        """.stripIndent()
    }
}
