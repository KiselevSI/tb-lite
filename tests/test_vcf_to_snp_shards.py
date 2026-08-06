import csv
import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "bin"))

from vcf_to_snp_shards import ALLELE_COLUMNS, SITE_COLUMNS, VCF_COLUMNS, extract_from_vcf

FIXTURE = REPO_ROOT / "tests" / "data" / "ERR4797591.mini.ann.vcf"

# Значения из продовой таблицы snp_sites (tb_database.dump). Совместимость по
# site_key обязательна: по нему новые данные схлопываются с существующими.
GOLDEN_SITE_KEYS = {
    1977: "7229c7508aa7fb285706a91c398ca4b22bdccee9",
    4013: "f03d041b126f9db3c71d4f62fb698395201f8ee5",
}


class ExtractFromVcfTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = extract_from_vcf(FIXTURE)

    def test_sample_id_from_vcf_header(self):
        self.assertEqual(self.result.sample_id, "ERR4797591")

    def test_site_key_matches_production(self):
        by_pos = {int(row[2]): row[0] for row in self.result.sites}
        for pos, expected in GOLDEN_SITE_KEYS.items():
            self.assertEqual(by_pos[pos], expected, f"site_key разошёлся на pos {pos}")

    def test_column_counts(self):
        self.assertTrue(all(len(row) == len(SITE_COLUMNS) for row in self.result.sites))
        self.assertTrue(all(len(row) == len(ALLELE_COLUMNS) for row in self.result.alleles))
        self.assertTrue(all(len(row) == len(VCF_COLUMNS) for row in self.result.vcf_rows))

    def test_keeps_every_variant_type(self):
        # 7 записей VCF, каждая с одной ALT -> 7 строк vcf_table
        self.assertEqual(len(self.result.vcf_rows), 7)
        self.assertEqual(
            {int(row[1]) for row in self.result.vcf_rows},
            {1977, 4013, 26747, 26957, 125830, 489935, 3133054},
        )

    def test_multi_ann_record_expands_sites_but_not_vcf_rows(self):
        # pos 489935 перекрывается двумя генами -> 2 site/allele, но 1 строка vcf_table
        sites_489935 = [row for row in self.result.sites if row[2] == "489935"]
        keys_489935 = {row[0] for row in sites_489935}
        alleles_489935 = [row for row in self.result.alleles if row[1] in keys_489935]
        vcf_489935 = [row for row in self.result.vcf_rows if row[1] == "489935"]
        self.assertEqual(len(sites_489935), 2)
        self.assertEqual(len(alleles_489935), 2)
        self.assertEqual(len(vcf_489935), 1)
        self.assertEqual(len(self.result.sites), 8)
        self.assertEqual(len(self.result.alleles), 8)

    def test_product_name_and_strand_always_empty(self):
        pn = SITE_COLUMNS.index("product_name")
        strand = SITE_COLUMNS.index("strand")
        self.assertTrue(all(row[pn] == "" for row in self.result.sites))
        self.assertTrue(all(row[strand] == "" for row in self.result.sites))

    def test_allele_metrics_match_production_semantics(self):
        af = ALLELE_COLUMNS.index("af")
        filt = ALLELE_COLUMNS.index("filter")
        dp = ALLELE_COLUMNS.index("dp")
        self.assertTrue(all(row[af] == "1" for row in self.result.alleles))
        self.assertTrue(all(row[filt] == "" for row in self.result.alleles))
        self.assertTrue(all(row[dp].isdigit() for row in self.result.alleles))

    def test_site_keys_unique_within_sample(self):
        keys = [row[0] for row in self.result.sites]
        self.assertEqual(len(keys), len(set(keys)))


if __name__ == "__main__":
    unittest.main()
