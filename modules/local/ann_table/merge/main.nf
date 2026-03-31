process MERGE {
  input:
    path(vcf)
    path(tbi)
  output:
    path "cohort.merge.vcf.gz", emit: vcf
  script:
    """
    bcftools merge -Oz -o cohort.merge.vcf.gz ${vcf}
    """
}