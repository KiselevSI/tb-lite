process RENAME_CHR {
    tag        "RENAME_CHR: $sample_name"
    

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