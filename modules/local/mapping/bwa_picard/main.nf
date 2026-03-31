process BWA_PICARD {
    tag "mapping: $sample_name"
    label 'big_mem'
    label 'multi_cpu'
    publishDir "${params.outdir}/mapped/${sample_name}", mode: params.mode

    input:
        tuple  val(sample_name), path(fastq_files)
        path   bwa_index
        path ref

    output:
    tuple val(sample_name), path("${sample_name}.dedup.bam"), path("${sample_name}.dedup.bam.bai"), emit: bam
    //path("*")

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