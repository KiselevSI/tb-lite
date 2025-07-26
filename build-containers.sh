singularity build --fakeroot containers/build_table.sif containers/def/build_table.def
singularity build --fakeroot containers/bwa-picard.sif containers/def/bwa-picard.def
singularity build --fakeroot containers/ismap.sif containers/def/ismap.def
singularity build --fakeroot containers/map.sif containers/def/map.def
singularity build --fakeroot containers/mapping.sif containers/def/mapping.def
singularity build --fakeroot containers/mosdpt.sif containers/def/mosdpt.def
singularity build --fakeroot containers/rd.sif containers/def/rd.def
singularity build --fakeroot containers/snpeff.sif containers/def/snpeff.def
singularity build --fakeroot containers/spotyping.sif containers/def/spotyping.def
singularity build --fakeroot containers/tb_mix.sif containers/def/tb_mix.def
singularity build --fakeroot containers/tb-lineage.sif containers/def/tb-lineage.def
singularity build --fakeroot containers/tblg.sif containers/def/tblg.def
singularity build --fakeroot containers/tbp.sif containers/def/tbp.def
singularity build --fakeroot containers/vcf2table.sif containers/def/vcf2table.def

singularity pull containers/fastqc_latest.sif docker://staphb/fastqc:0.12.1
singularity pull containers/bcftools.1.22.sif docker://staphb/bcftools:1.22
singularity pull containers/fastp_0.24.1.sif  docker://staphb/fastp:0.24.1
singularity pull containers/multiqc_v1.30.sif  docker://multiqc/multiqc:1.30