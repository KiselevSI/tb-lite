#!/usr/bin/env bash
# Склейка шардов SNP_DB_SHARDS в итоговые таблицы для импорта в TB Platform.
#
# На вход — каталог с sites.*.tsv.gz / alleles.*.tsv.gz / vcfrows.*.tsv.gz
# (без заголовков) плюс samples.*.txt и counts.*.tsv.
# На выход — три таблицы с заголовком, общий список образцов и манифест.
set -euo pipefail

SHARD_DIR=""
OUT_DIR=""
THREADS=4
SORT_MEM="2G"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --shard-dir) SHARD_DIR="$2"; shift 2 ;;
        --out-dir)   OUT_DIR="$2";   shift 2 ;;
        --threads)   THREADS="$2";   shift 2 ;;
        --sort-mem)  SORT_MEM="$2";  shift 2 ;;
        *) echo "Неизвестный аргумент: $1" >&2; exit 2 ;;
    esac
done

if [[ -z "$SHARD_DIR" || -z "$OUT_DIR" ]]; then
    echo "Использование: merge_snp_shards.sh --shard-dir DIR --out-dir DIR [--threads N] [--sort-mem 2G]" >&2
    exit 2
fi

mkdir -p "$OUT_DIR"
WORK_TMP="$(mktemp -d "${OUT_DIR}/merge_tmp.XXXXXX")"
trap 'rm -rf "$WORK_TMP"' EXIT

if command -v pigz >/dev/null 2>&1; then
    compress() { pigz -p "$THREADS"; }
else
    compress() { gzip; }
fi

SITES_HEADER=$'site_key\tchrom\tpos\tref\talt\tqual\teffect\timpact\tgene\tgene_id\tproduct_name\tstrand\tfeature_id\tbiotype\trank\thgvs_c\thgvs_p\tcds_pos\taa_pos'
ALLELES_HEADER=$'sample_id\tsite_key\tallele\tdp\taf\tqual\tfilter'
VCF_HEADER=$'sample_id\tpos\talt'

shopt -s nullglob
site_shards=( "$SHARD_DIR"/sites.*.tsv.gz )
allele_shards=( "$SHARD_DIR"/alleles.*.tsv.gz )
vcf_shards=( "$SHARD_DIR"/vcfrows.*.tsv.gz )
count_files=( "$SHARD_DIR"/counts.*.tsv )
sample_files=( "$SHARD_DIR"/samples.*.txt )
shopt -u nullglob

if (( ${#site_shards[@]} == 0 )); then
    echo "ERROR: в $SHARD_DIR нет шардов sites.*.tsv.gz" >&2
    exit 1
fi

# snp_sites: одна строка на site_key, побеждает максимальный QUAL — та же
# семантика, что была при первичной заливке прода. Поле 6 — qual, пустое
# значение считаем минимально возможным.
{
    printf '%s\n' "$SITES_HEADER"
    zcat "${site_shards[@]}" \
        | LC_ALL=C sort -t $'\t' -k1,1 -S "$SORT_MEM" --parallel="$THREADS" -T "$WORK_TMP" \
        | awk -F'\t' -v OFS='\t' '
            function q(v) { return (v == "" ? -1e308 : v + 0) }
            NR == 1 { key = $1; line = $0; best = q($6); next }
            $1 == key { if (q($6) > best) { line = $0; best = q($6) } ; next }
            { print line; key = $1; line = $0; best = q($6) }
            END { if (NR > 0) print line }
          '
} > "$WORK_TMP/snp_sites.tsv"

sites_rows=$(( $(wc -l < "$WORK_TMP/snp_sites.tsv") - 1 ))
compress < "$WORK_TMP/snp_sites.tsv" > "$OUT_DIR/snp_sites.tsv.gz"

{ printf '%s\n' "$ALLELES_HEADER"; zcat "${allele_shards[@]}"; } | compress > "$OUT_DIR/sample_snp_alleles.tsv.gz"
{ printf '%s\n' "$VCF_HEADER";     zcat "${vcf_shards[@]}";    } | compress > "$OUT_DIR/vcf_table.tsv.gz"

cat "${sample_files[@]}" | LC_ALL=C sort -u > "$OUT_DIR/samples.txt"

# Счётчики берём из counts.*.tsv шардов, чтобы не разжимать итоговые файлы ещё раз.
sum_metric() {
    local metric="$1"
    awk -F'\t' -v m="$metric" '$1 == m { total += $2 } END { print total + 0 }' "${count_files[@]}"
}

{
    printf 'shards\t%s\n'                  "${#site_shards[@]}"
    printf 'samples\t%s\n'                 "$(wc -l < "$OUT_DIR/samples.txt")"
    printf 'sites_rows_before_dedup\t%s\n' "$(sum_metric sites)"
    printf 'sites_rows\t%s\n'              "$sites_rows"
    printf 'alleles_rows\t%s\n'            "$(sum_metric alleles)"
    printf 'vcf_rows\t%s\n'                "$(sum_metric vcfrows)"
} > "$OUT_DIR/snp_db_manifest.tsv"

echo "OK: итоговые таблицы в $OUT_DIR"
