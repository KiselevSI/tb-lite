#!/usr/bin/env python3
"""Генерация psql-скрипта импорта SNP-таблиц в TB Platform."""
from __future__ import annotations

import argparse
from pathlib import Path

TEMPLATE = r"""-- Импорт SNP-таблиц TB Platform. Сгенерировано пайплайном tb-lite.
--
-- Запуск (файлы читает тот процесс, который выполняет psql, поэтому при
-- запуске psql в контейнере каталог с данными нужно смонтировать):
--     psql "$DATABASE_URL" -f import_snp.sql
--
-- Порядок важен: general.tsv должен быть импортирован раньше — на
-- sample_snp_alleles.sample_id висит FK на general."ID".
-- Повторный запуск безопасен на каждом шаге.

\set ON_ERROR_STOP on
\timing on

BEGIN;

\echo '=== staging ==='

CREATE TEMP TABLE stg_sites (
    site_key text, chrom text, pos integer, ref text, alt text, qual double precision,
    effect text, impact text, gene text, gene_id text, product_name text, strand text,
    feature_id text, biotype text, rank text, hgvs_c text, hgvs_p text,
    cds_pos text, aa_pos text
);
\copy stg_sites FROM PROGRAM 'zcat ''{data_dir}/snp_sites.tsv.gz''' WITH (FORMAT csv, DELIMITER E'\t', HEADER true, NULL '', QUOTE E'\b')

CREATE TEMP TABLE stg_alleles (
    sample_id text, site_key text, allele text, dp integer, af double precision,
    qual double precision, filter text
);
\copy stg_alleles FROM PROGRAM 'zcat ''{data_dir}/sample_snp_alleles.tsv.gz''' WITH (FORMAT csv, DELIMITER E'\t', HEADER true, NULL '', QUOTE E'\b')

CREATE TEMP TABLE stg_vcf (sample_id text, pos integer, alt text);
\copy stg_vcf FROM PROGRAM 'zcat ''{data_dir}/vcf_table.tsv.gz''' WITH (FORMAT csv, DELIMITER E'\t', HEADER true, NULL '', QUOTE E'\b')

\echo '=== snp_sites ==='

INSERT INTO snp_sites (
    site_key, chrom, pos, ref, alt, qual, effect, impact, gene, gene_id,
    product_name, strand, feature_id, biotype, rank, hgvs_c, hgvs_p, cds_pos, aa_pos
)
SELECT
    site_key, chrom, pos, ref, alt, qual, effect, impact, gene, gene_id,
    product_name, strand, feature_id, biotype, rank, hgvs_c, hgvs_p, cds_pos, aa_pos
FROM stg_sites
ON CONFLICT (site_key) DO NOTHING;

CREATE INDEX ON stg_alleles (site_key);
ANALYZE stg_alleles;

\echo '=== sample_snp_alleles ==='

INSERT INTO sample_snp_alleles (sample_id, site_id, allele, dp, af, qual, filter)
SELECT a.sample_id, s.site_id, a.allele, a.dp, a.af, a.qual, a.filter
FROM stg_alleles a
JOIN snp_sites s ON s.site_key = a.site_key
WHERE EXISTS (SELECT 1 FROM general g WHERE g."ID" = a.sample_id)
ON CONFLICT (sample_id, site_id) DO NOTHING;

\echo '=== vcf_table ==='

-- У vcf_table нет уникального ключа, поэтому повторный импорт — это
-- удаление строк входящих образцов и вставка заново.
CREATE TEMP TABLE stg_vcf_samples AS SELECT DISTINCT sample_id FROM stg_vcf;

DELETE FROM vcf_table v
USING stg_vcf_samples s
WHERE v.sample_id = s.sample_id;

INSERT INTO vcf_table (sample_id, pos, alt)
SELECT v.sample_id, v.pos, v.alt
FROM stg_vcf v
WHERE EXISTS (SELECT 1 FROM general g WHERE g."ID" = v.sample_id);

CREATE TEMP TABLE stg_new_samples AS
SELECT sample_id, ((row_number() OVER (ORDER BY sample_id)) - 1) / {chunk_size} AS chunk_no
FROM (
    SELECT DISTINCT a.sample_id
    FROM stg_alleles a
    WHERE EXISTS (SELECT 1 FROM general g WHERE g."ID" = a.sample_id)
) t;

CREATE INDEX ON stg_new_samples (chunk_no);
ANALYZE stg_new_samples;

COMMIT;

\echo '=== sample_snp_profiles ==='

-- Отдельно от загрузки данных: падение на профилях не должно откатывать
-- snp_sites/sample_snp_alleles/vcf_table. COMMIT внутри DO допустим, потому
-- что блок выполняется вне явной транзакции. Маски перебираются через
-- FOREACH по массиву, а не курсором: при открытом портале COMMIT запрещён.
DO $do$
DECLARE
    masks text[];
    mask_name text;
    chunk integer;
    max_chunk integer;
BEGIN
    SELECT array_agg(mask_key ORDER BY sort_order) INTO masks FROM snp_masks;
    IF masks IS NULL THEN
        RAISE NOTICE 'Таблица snp_masks пуста — профили не строятся. Загрузите маски (tb-platform/deploy/rebuild_snp_profiles_masked.sql).';
        RETURN;
    END IF;

    SELECT max(chunk_no) INTO max_chunk FROM stg_new_samples;
    IF max_chunk IS NULL THEN
        RAISE NOTICE 'Нет новых образцов — профили не строятся.';
        RETURN;
    END IF;

    FOREACH mask_name IN ARRAY masks LOOP
        DROP TABLE IF EXISTS tmp_canon;
        CREATE TEMP TABLE tmp_canon AS
        WITH snp AS (
            SELECT ss.site_id, ss.chrom, ss.pos, ss.ref, ss.alt
            FROM snp_sites ss
            WHERE char_length(ss.ref) = 1 AND char_length(ss.alt) = 1
              AND upper(ss.ref) IN ('A','C','G','T','N')
              AND upper(ss.alt) IN ('A','C','G','T','N')
              AND (
                    mask_name = 'none'
                    OR NOT EXISTS (
                        SELECT 1 FROM snp_mask_regions r
                        WHERE r.mask_key = mask_name AND r.pos_range @> ss.pos
                    )
              )
        )
        SELECT site_id, min(site_id) OVER (PARTITION BY chrom, pos, ref, alt) AS cid
        FROM snp;

        CREATE INDEX ON tmp_canon (site_id);
        ANALYZE tmp_canon;

        FOR chunk IN 0..max_chunk LOOP
            INSERT INTO sample_snp_profiles (mask_key, sample_id, site_ids, snp_count)
            SELECT mask_name,
                   a.sample_id,
                   array_agg(DISTINCT c.cid ORDER BY c.cid),
                   count(DISTINCT c.cid)
            FROM sample_snp_alleles a
            JOIN tmp_canon c ON c.site_id = a.site_id
            WHERE a.sample_id IN (
                SELECT sample_id FROM stg_new_samples WHERE chunk_no = chunk
            )
            GROUP BY a.sample_id
            ON CONFLICT (mask_key, sample_id) DO UPDATE
                SET site_ids = EXCLUDED.site_ids,
                    snp_count = EXCLUDED.snp_count,
                    updated_at = now();
            COMMIT;
        END LOOP;
    END LOOP;

    DROP TABLE IF EXISTS tmp_canon;
END
$do$;

\echo '=== ANALYZE ==='
ANALYZE snp_sites;
ANALYZE sample_snp_profiles;

\echo 'Готово. vcf_table и sample_snp_alleles стоит проанализировать отдельно:'
\echo '  ANALYZE vcf_table; ANALYZE sample_snp_alleles;'
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        required=True,
        help="Каталог с snp_sites.tsv.gz / sample_snp_alleles.tsv.gz / vcf_table.tsv.gz "
             "в той файловой системе, откуда будет запущен psql.",
    )
    parser.add_argument("-o", "--out", required=True, type=Path)
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=5000,
        help="Сколько образцов обрабатывать за одну транзакцию при построении профилей.",
    )
    args = parser.parse_args()

    data_dir = args.data_dir.rstrip("/")
    args.out.write_text(
        TEMPLATE.format(data_dir=data_dir, chunk_size=args.chunk_size),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
