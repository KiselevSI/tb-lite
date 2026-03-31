class WorkflowMain {

    public static void initialise(workflow, params, log) {

        // Print help message
        if (params.help) {
            log.info """
            ======================
            TB-Lite Pipeline v${workflow.manifest.version}
            ======================

            Usage:
              nextflow run main.nf --samples <samples.csv> [options]

            Required:
              --samples       CSV file with sample information (Sample,R1,R2,Layout)

            Options:
              --outdir        Output directory [default: ./results]
              --reference     Reference genome FASTA [default: assets/h37rv.fa]
              --skip_reports  Skip MultiQC and final reports [default: false]
              --mode          publishDir mode: copy or link [default: copy]
              --help          Show this help message

            Profiles:
              -profile local  Local execution (default)
              -profile k8s    Kubernetes execution
            """.stripIndent()
            System.exit(0)
        }

        // Validate required params
        if (!params.samples) {
            log.error "ERROR: --samples parameter is required. Use --help for usage information."
            System.exit(1)
        }

        // Log pipeline info
        log.info """
        ======================
        TB-Lite Pipeline v${workflow.manifest.version}
        ======================
        samples    : ${params.samples}
        outdir     : ${params.outdir}
        reference  : ${params.reference}
        profile    : ${workflow.profile}
        ======================
        """.stripIndent()
    }
}
