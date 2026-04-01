singularity build --fakeroot containers/ann_table.sif containers/def/ann_table.def
singularity build --fakeroot containers/build_table.sif containers/def/build_table.def
singularity build --fakeroot containers/tb_platform_tables.sif containers/def/tb_platform_tables.def
singularity build --fakeroot containers/rd.sif containers/def/rd.def
singularity build --fakeroot containers/spotyping.sif containers/def/spotyping.def
singularity build --fakeroot containers/tb_mix.sif containers/def/tb_mix.def
singularity build --fakeroot containers/tblg.sif containers/def/tblg.def
singularity pull containers/fastqc_latest.sif docker://staphb/fastqc:0.12.1
singularity pull containers/bcftools.1.22.sif docker://staphb/bcftools:1.22
singularity pull containers/fastp_0.24.1.sif  docker://staphb/fastp:0.24.1
singularity pull containers/multiqc_v1.30.sif  docker://multiqc/multiqc:v1.30
singularity pull containers/picard_3.4.0.sif  docker://broadinstitute/picard:3.4.0
