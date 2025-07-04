params.reads = "./data"
params.outdir = "./results"
params.fastqc_dir = "${params.outdir}/fastqc_reports"
params.threads = 3
params.reference    = "ref/h37rv.fa"
params.mode = 'link'
params.rd_path    = "scripts/rd.py"
params.rd_db    = "db/RD.bed"
params.IS6110    = "IS6110_db/is6110.fasta"
params.ref_gbk    = "IS6110_db/h37rv.gbk"

process run_fastqc {
    tag "fastqc: $sample_name"
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
    tag "fastp: $sample_name"
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
    tag "mapping: $sample_name"
    publishDir "${params.outdir}/mapped/${sample_name}", mode: params.mode

    input:
        tuple  val(sample_name), path(fastq_files)
        path   bwa_index
        path ref

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
    tag "call_variants: ${sample_name}"
    publishDir "${params.outdir}/vcf/${sample_name}", mode: params.mode

    input:
        tuple val(sample_name), path(bam), path(bam_idx)
        path ref
        

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
    tag "mosdepth: ${sample_name}"
    publishDir "${params.outdir}/stats/mosdepth/${sample_name}", mode: params.mode

    input:
        tuple val(sample_name), path(bam), path(bam_idx)        

    output:
        tuple path("${sample_name}.mosdepth.global.dist.txt"), path("${sample_name}.per-base.bed.gz.csi")

        tuple val(sample_name), path("${sample_name}.mosdepth.summary.txt"), emit: cov

        tuple val(sample_name), path("${sample_name}.per-base.bed.gz"),  emit: bed

        tuple val(sample_name), path("${sample_name}.median.txt"), emit: median



    script:
        """
        mosdepth --fast-mode $sample_name $bam

        sort -k2,2nr ${sample_name}.mosdepth.global.dist.txt | awk '\$3>=0.5 {print \$2; exit}' > ${sample_name}.median.txt

        """
}
process run_rd {
    tag "RD: ${sample_name}"
    publishDir "${params.outdir}/rd/${sample_name}", mode: params.mode

    input:
        tuple val(sample_name), path(bed)
        path rd
        path db

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

process run_spotyping {
    tag        "SpoTyping: $sample_name"
    publishDir "${params.outdir}/spotyping/$sample_name", mode: params.mode
    conda      "conda-envs/spotyping.yaml" 

    input:
        tuple val(sample_name), path(fastq_files)

    output:
        path("$sample_name.*"), emit: other
        path("$sample_name"), emit: code

    script:
        def reads = fastq_files instanceof List ? fastq_files : [fastq_files]

        reads = reads.join(" ")

        """
        SpoTyping.py $reads -o $sample_name
        """
}

process run_tblg {
    tag        "tblg: $sample_name"
    publishDir "${params.outdir}/lineage", mode: params.mode
    conda      "conda-envs/tblg.yaml" 

    input:
        tuple val(sample_name), path(vcf), path(vcf_csi)

    output:
        path("${sample_name}.lg.tsv")

    script:
        
        """
        tblg $vcf -o ${sample_name}.lg.tsv
        """
}

process run_is6110 {
    tag        "is6110: $sample_name"
    publishDir "${params.outdir}/is6110/paired", mode: params.mode
    conda      "conda-envs/ismapper.yaml" 

    input:
        tuple val(sample_name), path(read1), path(read2)
        each IS6110
        each ref_gbk

    output:
        path("$sample_name/*")

    script:

        """
        ismap --reads $read1 $read2 --queries $IS6110 --reference $ref_gbk --bam
        """
}

process run_rename_chromosome {
    tag        "drug_resist: $sample_name"
    

    input:
        tuple val(sample_name), path(vcf), path(vcf_csi)
        path chromosome_name
        

    output:
        tuple val(sample_name), path("${sample_name}.renamed_chromosome.vcf.gz")

    script:

        """
        bcftools annotate --rename-chrs $chromosome_name $vcf -O z -o ${sample_name}.renamed_chromosome.vcf.gz
        """
}

process run_annotate_vcf{
    tag        "drug_resist: $sample_name"
    conda      "conda-envs/dr.yaml" 

    input:
        tuple val(sample_name), path(vcf_renamed)        

    output:
        tuple val(sample_name), path("${sample_name}.annotated.vcf.gz")

    script:

        """
        snpEff ann -v Mycobacterium_tuberculosis_h37rv $vcf_renamed | bgzip -c > ${sample_name}.annotated.vcf.gz
        """
}

process run_drug_resist {
    tag        "drug_resist: $sample_name"
    publishDir "${params.outdir}/drug_resist/$sample_name", mode: params.mode
    conda      "conda-envs/dr.yaml" 

    input:
        tuple val(sample_name), path(vcf_annotated)
        path tb_resistance
        path db_drug_resist
        

    output:
        path("*")

    script:

        """
        python3 $tb_resistance -i $vcf_annotated -o ${sample_name}.drug_resist.csv -d
        """
}


workflow {

    reads = Channel.fromPath("${params.reads}/*.fastq", checkIfExists: true)
        .map { [it.baseName.replaceFirst(/_R?[12].*/, ''), it] }
        .groupTuple()

    

    IS6110 = Channel.fromPath(params.IS6110)
    ref_gbk = Channel.fromPath(params.ref_gbk)


        
    

    trimmed = run_fastp(reads).trimmed_reads

    bwa_idx = Channel.fromPath("${params.reference}.*")          // ref/h37rv.fa.*
          .filter { !it.name.endsWith('.fa') }               // вдруг попала копия .fa
          .collect()

    ref = Channel.value(file("$params.reference"))

    bam_all = run_mapping(trimmed, bwa_idx, ref)

    cov_all = run_mosdepth(bam_all)
    
    

    // ----------------------------------------------------------------
    // 4) Определяем «хорошие» образцы (медиана > 0)
    // ----------------------------------------------------------------
    cov_all.median
        .map    { id, f -> tuple(id, ((f.text.trim()) ?: '0') as Integer) }
        .filter { id, cov -> cov > 5 }
        .map    { id, cov -> id }  // берём только первый элемент
        .set    { GOOD_SAMPLES }         // дублируемый канал


    // «плохие» ─ медиана ≤ 0
    cov_all.median
        .map    { id, f -> tuple(id, ((f.text.trim()) ?: '0') as Integer) }
        .filter { id, cov -> cov <= 5 }
        .map    { id, cov -> id }
        .set    { BAD_SAMPLES }       // новый канал


    // results/bad_samples.txt с перезаписью при повторных запусках
    BAD_SAMPLES
        .collectFile(
            name:      'bad_samples_coverage_less_than_5.txt',   // имя результирующего файла
            storeDir:  params.outdir,       // куда положить
            newLine:   true,                // добавлять \n между элементами
            sort:      true               // (опционно) отсортировать
            )






    // ----------------------------------------------------------------
    // 5) Фильтруем trimmed‑reads, BAM и BED
    // ----------------------------------------------------------------
    trimmed_good = trimmed
        .join( GOOD_SAMPLES.map { tuple(it) } )
        .map  { id, fq_files -> tuple(id, fq_files) }

    bam_good = bam_all
        .join( GOOD_SAMPLES.map { tuple(it) } )
        .map  { id, bam_file, bai_file -> tuple(id, bam_file, bai_file) }

    bed_good = cov_all.bed
        .join( GOOD_SAMPLES.map { tuple(it) } )
        .map  { id, bed_file -> tuple(id, bed_file) }


    run_fastqc(trimmed_good)

    run_spotyping(trimmed_good)

    //single_reads =  reads.filter { id, files -> files.size() == 1 }
    //    .map    { id, files -> tuple(id, files[0]) }.groupTuple()
    
    paired_reads = reads.filter { _id, files -> files.size() == 2 }
        .map    { id, files ->
            // упорядочим, чтобы сначала был R1
            def (r1, r2) = files.sort { it.name }
            tuple(id, r1, r2)
        }.join( GOOD_SAMPLES.map { tuple(it) } )

    run_is6110(paired_reads, IS6110, ref_gbk)


    rd_path = Channel.value(file(params.rd_path))
    rd_db = Channel.value(file(params.rd_db))

    

    run_rd(cov_all.bed, rd_path, rd_db)

    vcf = run_call_variants(bam_good, ref)

    db_drug_resist = Channel.value(file("db_drug_resist"))
    chr_name = Channel.value(file("scripts/chr.txt"))
    script_dr_path = Channel.value(file("scripts/tb_resistance.py"))

    vcf_renamed = run_rename_chromosome(vcf, chr_name)
    vcf_annotated = run_annotate_vcf(vcf_renamed)

    run_drug_resist(vcf_annotated, script_dr_path, db_drug_resist)

    run_tblg(vcf)
}