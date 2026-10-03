import csv
import gzip
import subprocess
import sys
import tempfile
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


def _read_gz_rows(path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return [row for row in csv.reader(handle, delimiter="\t") if row]


def _make_variant_vcf(target: Path, sample: str, qual_factor: float) -> Path:
    """Копия фикстуры с другим именем образца и изменённым QUAL."""
    lines = []
    for line in FIXTURE.read_text(encoding="utf-8").splitlines():
        if line.startswith("#CHROM"):
            lines.append(line.replace("ERR4797591", sample))
        elif line.startswith("#"):
            lines.append(line)
        else:
            parts = line.split("\t")
            parts[5] = f"{float(parts[5]) * qual_factor:.2f}"
            lines.append("\t".join(parts))
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


class ShardCliTest(unittest.TestCase):
    """Два образца с одинаковыми сайтами, но разным QUAL: проверяем дедупликацию."""

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        tmp = Path(cls._tmp.name)

        second = _make_variant_vcf(tmp / "SAMPLE2.annotated.ann.vcf", "SAMPLE2", 0.5)

        file_list = tmp / "vcfs.txt"
        file_list.write_text(f"{FIXTURE}\n{second}\n", encoding="utf-8")

        cls.outdir = tmp / "out"
        subprocess.run(
            [sys.executable, str(REPO_ROOT / "bin" / "vcf_to_snp_shards.py"),
             "--file-list", str(file_list), "--tag", "batch_1",
             "-o", str(cls.outdir), "-j", "2"],
            check=True,
        )

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_writes_expected_files(self):
        for name in ("sites.batch_1.tsv.gz", "alleles.batch_1.tsv.gz",
                     "vcfrows.batch_1.tsv.gz", "samples.batch_1.txt",
                     "counts.batch_1.tsv"):
            self.assertTrue((self.outdir / name).exists(), name)

    def test_shards_have_no_header(self):
        first = _read_gz_rows(self.outdir / "alleles.batch_1.tsv.gz")[0]
        self.assertNotEqual(first[0], "sample_id")

    def test_row_counts(self):
        # 8 уникальных site_key на два образца, по 8 аллелей и 7 vcf-строк на образец
        self.assertEqual(len(_read_gz_rows(self.outdir / "sites.batch_1.tsv.gz")), 8)
        self.assertEqual(len(_read_gz_rows(self.outdir / "alleles.batch_1.tsv.gz")), 16)
        self.assertEqual(len(_read_gz_rows(self.outdir / "vcfrows.batch_1.tsv.gz")), 14)

    def test_sites_keep_max_qual(self):
        rows = _read_gz_rows(self.outdir / "sites.batch_1.tsv.gz")
        by_pos = {int(row[2]): float(row[5]) for row in rows}
        # У ERR4797591 QUAL вдвое выше, значит должен победить он
        self.assertAlmostEqual(by_pos[1977], 1643.47, places=2)

    def test_samples_file(self):
        samples = (self.outdir / "samples.batch_1.txt").read_text().split()
        self.assertEqual(sorted(samples), ["ERR4797591", "SAMPLE2"])

    def test_counts_file(self):
        counts = dict(
            line.split("\t")
            for line in (self.outdir / "counts.batch_1.tsv").read_text().splitlines()
        )
        self.assertEqual(counts["sites"], "8")
        self.assertEqual(counts["alleles"], "16")
        self.assertEqual(counts["vcfrows"], "14")
        self.assertEqual(counts["samples"], "2")


if __name__ == "__main__":
    unittest.main()
