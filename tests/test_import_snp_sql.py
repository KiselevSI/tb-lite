import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
GENERATOR = REPO_ROOT / "bin" / "write_import_snp_sql.py"


class WriteImportSnpSqlTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        out = Path(cls._tmp.name) / "import_snp.sql"
        subprocess.run(
            [sys.executable, str(GENERATOR), "--data-dir", "/data/snp",
             "-o", str(out), "--chunk-size", "5000"],
            check=True,
        )
        cls.sql = out.read_text(encoding="utf-8")

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_stops_on_error(self):
        self.assertIn("\\set ON_ERROR_STOP on", self.sql)

    def test_copies_all_three_files(self):
        self.assertIn("/data/snp/snp_sites.tsv.gz", self.sql)
        self.assertIn("/data/snp/sample_snp_alleles.tsv.gz", self.sql)
        self.assertIn("/data/snp/vcf_table.tsv.gz", self.sql)

    def test_idempotent_inserts(self):
        self.assertIn("ON CONFLICT (site_key) DO NOTHING", self.sql)
        self.assertIn("ON CONFLICT (sample_id, site_id) DO NOTHING", self.sql)
        self.assertIn("ON CONFLICT (mask_key, sample_id) DO UPDATE", self.sql)

    def test_vcf_table_is_delete_then_insert(self):
        self.assertIn("DELETE FROM vcf_table v", self.sql)

    def test_only_samples_present_in_general(self):
        self.assertIn('FROM general g WHERE g."ID"', self.sql)

    def test_profiles_use_canonical_site_id(self):
        self.assertIn("min(site_id) OVER (PARTITION BY chrom, pos, ref, alt)", self.sql)

    def test_chunk_size_applied(self):
        self.assertIn("/ 5000", self.sql)

    def test_data_load_is_one_transaction_and_profiles_are_not(self):
        self.assertLess(self.sql.index("BEGIN;"), self.sql.index("COMMIT;"))
        self.assertGreater(self.sql.index("DO $do$"), self.sql.index("COMMIT;"))


if __name__ == "__main__":
    unittest.main()
