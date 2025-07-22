
params.outdir = "./results"
params.reference    = "ref/h37rv.fa"
params.mode = 'link'
params.rd_db    = "db/RD.bed"
params.is6110    = "IS6110_db/is6110.fasta"
params.ref_gbk    = "IS6110_db/h37rv.gbk"
params.chr_name = "scripts/chr.txt"

nextflow.enable.dsl = 2 

params.multiqc = "multiqc_config.yaml"

params.min_align_pct  = 80
params.min_mean_cov   = 10
params.min_median   = 10

params.samples = "sample.csv"

workflow {


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

    fastqc_reports = run_fastqc(reads)
    
    trimmed = run_fastp(reads).trimmed_reads

    bwa_idx = Channel.fromPath("${params.reference}.*").filter { !it.name.endsWith('.fa') }.collect()

    ref = Channel.value(file("$params.reference"))

    bam_all = run_mapping(trimmed, bwa_idx, ref).bam

    levels = Channel.value(file("scripts/levels.tsv"))

    run_tb_mix(bam_all, ref, levels)


    map_stats = run_map_stats(bam_all, ref)

    wgs_metrics = map_stats.wgs
    align_metrics = map_stats.align
    samtools_stats = run_samtools_stats(bam_all)
    // ------------------ 1) процент выравненных ------------------
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



    // объединяем по ключу-идентификатору
    // good_samples = align_pct.join(mean_cov).join(median_cov)
    //     .filter { _id, pct, cov -> pct >= params.min_align_pct &&
    //                             cov >= params.min_mean_cov }
    //     .map { id, _p, _c -> id }

    good_samples = median_cov.filter{_id, median -> median >= params.min_median}.map {id, m -> id}
    

    trimmed_good  = trimmed.join(good_samples)

    bam_good      = bam_all.join(good_samples)

    mosdepth_coverage = run_mosdepth(bam_good)
    

    run_spotyping(trimmed_good)

    paired_reads = trimmed_good.filter { _id, files -> files.size() == 2 }
            .join( good_samples)


    is6110 = Channel.value(file(params.is6110))
    ref_gbk = Channel.value(file(params.ref_gbk))

    run_is6110(paired_reads, is6110, ref_gbk)


    rd_db = Channel.value(file(params.rd_db))

    

    run_rd(mosdepth_coverage.bed, rd_db)

    vcf = run_call_variants(bam_good, ref)

    run_tb_profiler_dr(vcf)

    chr_name = Channel.value(file(params.chr_name))

    vcf_renamed = run_rename_chromosome(vcf, chr_name)
    vcf_annotated = run_annotate_vcf(vcf_renamed)

    //run_drug_resist(vcf_annotated)

    run_tblg(vcf)

    cfg = Channel.fromPath(params.multiqc)

    bcftools_stats = run_bcftools_stats(vcf)

    multiqc_files = wgs_metrics.mix(align_metrics).mix(fastqc_reports).mix(bcftools_stats).mix(samtools_stats.stats).mix(samtools_stats.flagstat)

    run_multiqc(multiqc_files.collect().ifEmpty([]), cfg)

}




process run_fastqc {
    tag "fastqc: $sample_name"
    publishDir("${params.outdir}/fastqc/", mode: params.mode)
    
    input:
    tuple val(sample_name), path(fastq_files)
    
    output:
    path "*"
    
    script:

    def files = fastq_files instanceof List ? fastq_files : [fastq_files]

    if (files.size() == 2){

        def read1 = files[0]
        def read2 = files[1]
        """
        mkdir -p $sample_name
        fastqc -t ${task.cpus} -o $sample_name $read1 $read2
        """
    }else{
        def read1 = files[0]
        """
        mkdir -p $sample_name
        fastqc -t ${task.cpus} -o $sample_name $read1
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
    tuple val(sample_name), path("*_trimmed*fastq.gz"), emit: trimmed_reads    
    
    script:
    // Приводим fastq_files к списку, если это не список
    def files = fastq_files instanceof List ? fastq_files : [fastq_files]
    if (files.size() == 2) {
        def read1 = files[0]
        def read2 = files[1]
        """
        fastp -i $read1 -I $read2 -o ${sample_name}_trimmed_1.fastq.gz -O ${sample_name}_trimmed_2.fastq.gz --detect_adapter_for_pe --trim_poly_g
        """
    } else if (files.size() == 1) {
        def read1 = files[0]
        """
        fastp -i $read1 -o ${sample_name}_trimmed_1.fastq.gz --trim_poly_g
        """
    }
}

process run_mapping2 {
    tag "mapping: $sample_name"
    publishDir "${params.outdir}/mapped/${sample_name}", mode: params.mode

    input:
        tuple  val(sample_name), path(fastq_files)
        path   bwa_index
        path ref

    output:
    tuple val(sample_name), path("${sample_name}.dedup.bam"), path("${sample_name}.dedup.bam.bai"), emit: bam
    path("*")

    script:

        def files = fastq_files instanceof List ? fastq_files : [fastq_files]
            if (files.size() == 2) {
                def read1 = files[0]
                def read2 = files[1]
                """
                   bwa mem -t ${task.cpus} \
                     -R "@RG\\tID:${sample_name}\\tSM:${sample_name}\\tPL:ILLUMINA" \
                     ${ref} ${read1} ${read2} \
                    | samblaster \
                        --removeDups \
                        --addMateTags \
                    | samtools sort -@ ${task.cpus} -o ${sample_name}.dedup.bam - \
                    && samtools index ${sample_name}.dedup.bam
                """
                
            } else if (files.size() == 1) {
                def read1 = files[0]
                """
                   bwa mem -t ${task.cpus} \
                     -R "@RG\\tID:${sample_name}\\tSM:${sample_name}\\tPL:ILLUMINA" \
                     ${ref} ${read1}  \
                    | samblaster \
                        --removeDups \
                        --addMateTags \
                    | samtools sort -@ ${task.cpus} -o ${sample_name}.dedup.bam - \
                    && samtools index ${sample_name}.dedup.bam
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
    tuple val(sample_name), path("${sample_name}.dedup.bam"), path("${sample_name}.dedup.bam.bai"), emit: bam
    path("*")

    script:

        def files = fastq_files instanceof List ? fastq_files : [fastq_files]
            if (files.size() == 2) {
                def read1 = files[0]
                def read2 = files[1]
                """
                   bwa mem -M -t ${task.cpus} \
                     -R "@RG\\tID:${sample_name}\\tSM:${sample_name}\\tPL:ILLUMINA" \
                     ${ref} ${read1} ${read2} \
                    | samtools view -Sb - \
                    | samtools sort -@ ${task.cpus} -o ${sample_name}.sorted.bam

                    java -jar /usr/local/bin/picard.jar MarkDuplicates \
                        INPUT=${sample_name}.sorted.bam\
                        OUTPUT=${sample_name}.dedup.bam \
                        METRICS_FILE=dup_metrics.txt \
                        ASSUME_SORTED=true \
                        VALIDATION_STRINGENCY=SILENT



                    samtools index ${sample_name}.dedup.bam

                    rm ${sample_name}.sorted.bam
                """
                
            } else if (files.size() == 1) {
                def read1 = files[0]
                """
                   bwa mem -M -t ${task.cpus} \
                     -R "@RG\\tID:${sample_name}\\tSM:${sample_name}\\tPL:ILLUMINA" \
                     ${ref} ${read1}  \
                    | samtools view -Sb - \
                    | samtools sort -@ ${task.cpus} -o ${sample_name}.sorted.bam

                    java -jar /usr/local/bin/picard.jar MarkDuplicates \
                        INPUT=${sample_name}.sorted.bam\
                        OUTPUT=${sample_name}.dedup.bam \
                        METRICS_FILE=dup_metrics.txt \
                        ASSUME_SORTED=true \
                        VALIDATION_STRINGENCY=SILENT


                    samtools index ${sample_name}.dedup.bam

                    rm ${sample_name}.sorted.bam
                """
            }

}

process run_tb_mix {
    tag "run_lineage: ${sample_name}"
    publishDir "${params.outdir}/tb-mix", mode: params.mode

    input:
        tuple val(sample_name), path(bam), path(bai)
        path ref
        path levels

    output:
        path("${sample_name}.mix.tsv")

    script:
        """
        tb_mix.py -i $bam -r $ref -l $levels --mq 30 --bq 20 -f 0.05 -p $sample_name -o ${sample_name}.mix.tsv

        """
}



process run_map_stats {
    tag "map_stats: ${sample_name}"
    publishDir "${params.outdir}/stats/picard/${sample_name}", mode: params.mode

    input:
        tuple val(sample_name), path(bam), path(bai)
        path ref

    output:
        path("${sample_name}.wgs_metrics.txt"), emit: wgs
        path("${sample_name}.alignment_metrics.txt"), emit: align
        tuple val(sample_name), path("${sample_name}.mean_coverage.txt"), emit: mean
        tuple val(sample_name), path("${sample_name}.median_coverage.txt"), emit: median

    script:
        """
        java -jar /usr/picard/picard.jar \
        CollectWgsMetrics \
        I=$bam \
        O=${sample_name}.wgs_metrics.txt \
        R=$ref COUNT_UNPAIRED=true

        java -jar /usr/picard/picard.jar \
          CollectAlignmentSummaryMetrics \
          R=$ref \
          I=$bam \
          O=${sample_name}.alignment_metrics.txt
        grep -v "^#" ${sample_name}.wgs_metrics.txt | awk 'NR==3 {print \$2}'> ${sample_name}.mean_coverage.txt
        grep -v "^#" ${sample_name}.wgs_metrics.txt | awk 'NR==3 {print \$4}' > ${sample_name}.median_coverage.txt
        """

}

process run_samtools_stats{
    tag "run_samtools_stats: ${sample_name}"
    publishDir "${params.outdir}/stats/samtools/${sample_name}", mode: params.mode

    input:
        tuple val(sample_name), path(bam), path(bai)

    output:
        path("${sample_name}.samtools.txt"), emit: stats
        path("${sample_name}.flagstat.txt"),  emit: flagstat
        tuple val(sample_name), path("${sample_name}.rmp.txt"), emit: rmp

    script:
        """
        samtools stats ${bam} > ${sample_name}.samtools.txt
        samtools flagstat ${bam} > ${sample_name}.flagstat.txt

        perc=\$(grep -oP '\\(\\K[0-9]+\\.[0-9]+(?=%)' \
            ${sample_name}.flagstat.txt | head -n1)
        if [ -z \"\$perc\" ]; then perc="0.0"; fi
        echo \"\$perc\" > ${sample_name}.rmp.txt

        """
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
        bcftools mpileup --threads ${task.cpus} \
            --min-MQ 10 --ignore-overlaps --max-depth 3000 \
            -f ${ref} ${bam} -Ou | \
        bcftools call --threads ${task.cpus} --multiallelic-caller \
            --ploidy 1 --variants-only -Ou - | \
        bcftools view --threads ${task.cpus} \
            --include 'DP>10' -Oz -o ${sample_name}.vcf.gz
        bcftools index ${sample_name}.vcf.gz
        """
}

process run_bcftools_stats{
    tag "run_bcftools_stats: ${sample_name}"
    publishDir "${params.outdir}/stats/bcftools", mode: params.mode

    input:
        tuple val(sample_name), path(vcf), path(vcf_csi)

    output:
        path("${sample_name}.bcftools.txt")

    script:
        """
        bcftools stats $vcf > ${sample_name}.bcftools.txt
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
        mosdepth -t ${task.cpus} $sample_name $bam

        sort -k2,2nr ${sample_name}.mosdepth.global.dist.txt | awk '\$3>=0.5 {print \$2; exit}' > ${sample_name}.median.txt

        """
}
process run_rd {
    tag "RD: ${sample_name}"
    publishDir "${params.outdir}/rd/${sample_name}", mode: params.mode

    input:
        tuple val(sample_name), path(bed)
        path db

    output:
        tuple val(sample_name), path("${sample_name}.novel_rd.tsv"),
        path("${sample_name}.known_rd.tsv")

    script:
        """
        rd.py $bed \
        -k $db \
        -n ${sample_name}.novel_rd.tsv \
        -o ${sample_name}.known_rd.tsv
        """
}

process run_spotyping {
    tag        "SpoTyping: $sample_name"
    publishDir "${params.outdir}/spotyping/$sample_name", mode: params.mode

    input:
        tuple val(sample_name), path(fastq_files)

    output:
        path("$sample_name.*"), emit: other
        path("$sample_name"), emit: code

    script:
        def reads = fastq_files instanceof List ? fastq_files : [fastq_files]

        reads = reads.join(" ")

        """
        SpoTyping.py $reads --noQuery -o $sample_name
        awk -F '\\t' -v OFS='\\t' -v S='${sample_name}' \\
            'NR==1 {print S, \$2, \$3}' \\
            ${sample_name} > ${sample_name}.tsv
        """
}

process run_tblg {
    tag        "tblg: $sample_name"
    publishDir "${params.outdir}/lineage", mode: params.mode

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

    input:
        tuple val(sample_name), path(fastq_files)
        path is6110
        path ref_gbk

    output:
        path("${sample_name}/*")

    script:

        def reads = fastq_files instanceof List ? fastq_files : [fastq_files]
        reads = reads.join(" ")

        """
        ismap --reads $reads --queries $is6110 --reference $ref_gbk --bam --t ${task.cpus} --output_dir $sample_name
        """
}


process run_tb_profiler_dr {
    tag "run_tb_profiler_dr: $sample_name"
    publishDir "${params.outdir}/tb-profiler/drug-resist/$sample_name", mode: params.mode

    input:
        tuple val(sample_name), path(vcf), path(vcf_csi)

    output:
        path("*")

    script:

        """
        tb-profiler profile --vcf $vcf -p $sample_name --csv --txt
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

    input:
        tuple val(sample_name), path(vcf_renamed)        

    output:
        tuple val(sample_name), path("${sample_name}.annotated.vcf.gz")

    script:

        """
        java -jar /snpEff.jar ann -v Mycobacterium_tuberculosis_h37rv $vcf_renamed | bgzip -c > ${sample_name}.annotated.vcf.gz
        """
}

process run_drug_resist {
    tag        "drug_resist: $sample_name"
    publishDir "${params.outdir}/drug_resist/$sample_name", mode: params.mode

    input:
        tuple val(sample_name), path(vcf_annotated)
   

    output:
        path("*")

    script:

        """
        tb_resistance.py -i $vcf_annotated -o ${sample_name}.drug_resist.csv -d
        """
}








process run_multiqc {
    tag "multiqc"
    publishDir "${params.outdir}/multiqc", mode: params.mode

    input:
    path reports
    path cfg

    output:
    path "*"

    script:
    """
    multiqc --config $cfg --title 'TB-Lite QC' --filename tb_multiqc_report .
    """
}


