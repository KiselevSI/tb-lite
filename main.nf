params.reads = "./data"
params.outdir = "./results"
params.fastqc_dir = "${params.outdir}/fastqc_reports"
params.threads = 10
params.reference    = "ref/h37rv.fa"
params.mode = 'link'
params.rd_path    = "scripts/rd.py"
params.rd_db    = "/home/zerg/git/tb-lite/db/RD.bed"

process run_fastqc {
    tag "Sample: $sample_name"
    publishDir("${params.fastqc_dir}/$sample_name", mode: params.mode)
    
    input:
    tuple val(sample_name), path(fastq_files)
    
    output:
    path "*.{html,zip}"
    
    script:

    def files = fastq_files instanceof List ? fastq_files : [fastq_files]

    if (files.size() == 2){

        def read1 = files[0]
        def read2 = files[1]
        """
        fastqc -t ${params.threads} -o ./ $read1 $read2
        """
    }else{
        def read1 = files[0]
        """
        fastqc -t ${params.threads} -o ./ $read1
        """
    }


}

process run_fastp {
    tag "Sample: $sample_name"
    publishDir("${params.outdir}/fastp/$sample_name", mode: params.mode)
    
    input:
    tuple val(sample_name), path(fastq_files)
    
    output:
    path "*.{html,json}", emit: report_trimm
    tuple val(sample_name), path("*.trimmed.fastq"), emit: trimmed_reads    
    
    script:
    // Приводим fastq_files к списку, если это не список
    def files = fastq_files instanceof List ? fastq_files : [fastq_files]
    if (files.size() == 2) {
        def read1 = files[0]
        def read2 = files[1]
        """
        fastp -i $read1 -I $read2 -o ${sample_name}_1.trimmed.fastq -O ${sample_name}_2.trimmed.fastq --detect_adapter_for_pe --trim_poly_g
        """
    } else if (files.size() == 1) {
        def read1 = files[0]
        """
        fastp -i $read1 -o ${sample_name}_1.trimmed.fastq --detect_adapter_for_pe --trim_poly_g
        """
    }
}

process run_mapping {
    tag "Sample: $sample_name"
    publishDir "${params.outdir}/mapped/${sample_name}", mode: params.mode

    input:
        tuple  val(sample_name), path(fastq_files)
        path   bwa_index
        each ref

    output:
        tuple val(sample_name), path("${sample_name}.bam"), path("${sample_name}.bam.bai")

    script:
        

        def files = fastq_files instanceof List ? fastq_files : [fastq_files]
        //def ref = bwa_index instanceof List ? bwa_index : [bwa_index]
        //ref = ref[0]

        if (files.size() == 2) {
            def read1 = files[0]
            def read2 = files[1]
            """
            bwa mem -t ${params.threads} ${ref} ${read1} ${read2} \
              | samtools view -bS | samtools sort -o ${sample_name}.bam

            samtools index ${sample_name}.bam
            """
        } else if (files.size() == 1) {
            def read1 = files[0]
            """
            bwa mem -t ${params.threads} ${ref} ${read1} \
              | samtools view -bS | samtools sort -o  ${sample_name}.bam

            samtools index ${sample_name}.bam
            """
        }
}
process run_call_variants {
    tag "Sample: ${sample_name}"
    publishDir "${params.outdir}/vcf/${sample_name}", mode: params.mode

    input:
        tuple val(sample_name), path(bam), path(bam_idx)
        each ref
        

    output:
        tuple val(sample_name), path("${sample_name}.vcf.gz"), path("${sample_name}.vcf.gz.csi")

    script:
        """
        bcftools mpileup --threads ${params.threads} \
            --min-MQ 30 --ignore-overlaps --max-depth 3000 \
            -f ${ref} ${bam} -Ou | \
        bcftools call --threads ${params.threads} --multiallelic-caller \
            --ploidy 1 --variants-only -Ou - | \
        bcftools view --threads ${params.threads} \
            --include 'QUAL>20 && DP>10' -Oz -o ${sample_name}.vcf.gz
        bcftools index ${sample_name}.vcf.gz
        """
}

process run_mosdepth {
    tag "Sample: ${sample_name}"
    publishDir "${params.outdir}/stats/mosdepth/${sample_name}", mode: params.mode

    input:
        tuple val(sample_name), path(bam), path(bam_idx)        

    output:
        tuple path("${sample_name}.mosdepth.global.dist.txt"), path("${sample_name}.mosdepth.summary.txt"), path("${sample_name}.per-base.bed.gz.csi")

        tuple val(sample_name), path("${sample_name}.per-base.bed.gz"),  emit: bed

    script:
        """
        mosdepth --fast-mode $sample_name $bam
        """
}
process run_rd {
    tag "Sample: ${sample_name}"
    publishDir "${params.outdir}/rd/${sample_name}", mode: params.mode

    input:
        tuple val(sample_name), path(bed)
        each rd
        each db

    output:
        tuple val(sample_name), path("${sample_name}.novel_rd.tsv"),
        path("${sample_name}.known_rd.tsv")

    script:
        """
        python3 $rd $bed \
        -k $db \
        -n ${sample_name}.novel_rd.tsv \
        -o ${sample_name}.known_rd.tsv
        """
}
workflow {

    reads = Channel.fromPath("${params.reads}/*.fastq", checkIfExists: true)
        .map { [it.baseName.replaceFirst(/_R?[12].*/, ''), it] }
        .groupTuple()
        
    run_fastqc(reads)
    trimmed = run_fastp(reads).trimmed_reads

    bwa_idx = Channel.fromPath("${params.reference}*").collect()

    ref = Channel.fromPath("$params.reference")

    bam = run_mapping(trimmed, bwa_idx, ref)

    stat_mosdp = run_mosdepth(bam).bed

    rd_path = Channel.fromPath(params.rd_path)
    rd_db = Channel.fromPath(params.rd_db)

    run_rd(stat_mosdp, rd_path, rd_db)

    //ref = Channel.fromPath("$params.reference")

    run_call_variants(bam, ref)
}