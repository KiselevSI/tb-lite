

include { MAKE_TABLE } from '../../modules/local/ann_table/make_table/main'
include { MERGE } from '../../modules/local/ann_table/merge/main'
include { SNPEFF_ANN } from '../../modules/local/ann_table/snpeff_ann/main'
include { POST_PROCESS_TABLE } from '../../modules/local/ann_table/post_process_table/main'

workflow ANN_TABLE {
    take:
        vcfs

    main:
        chr_name = Channel.value(file(params.chr_name))
        other_count = vcfs.map { _sample_name, vcf, tbi, count ->
            tuple( vcf, tbi, count.text.trim() as Float )
        }

        filtered = other_count.filter { _vcf, _tbi, count -> count < 5000 }

        vcf = filtered.map { it[0]}.collect()
        tbi = filtered.map { it[1]}.collect()

        merged = MERGE(vcf, tbi)
        ann = SNPEFF_ANN(merged, chr_name)
        table = MAKE_TABLE(ann)
        feature_table = Channel.fromPath(params.h37rv_feature_table)
        final_table = POST_PROCESS_TABLE(table, feature_table)

}
