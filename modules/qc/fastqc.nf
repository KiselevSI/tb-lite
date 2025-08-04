process FASTQC {
    tag "fastqc: $sample_name"
    label 'medium_mem'
    label 'multi_cpu'
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