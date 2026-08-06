import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
import uuid
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
GENERATOR = REPO_ROOT / "bin" / "write_import_snp_sql.py"

# Интеграционный тест поднимает СВОЙ одноразовый Postgres и удаляет его в
# teardown. К рабочим базам он не подключается: у tb_platform_db на дев-машине
# лежит полная копия боевых данных, и любой незаскоупленный DELETE там
# необратим.
PG_IMAGE = "postgres:17"


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


# Схема, которой достаточно для проверки import_snp.sql. Колонки и ограничения
# скопированы с продовой схемы TB Platform (backend/db.py + миграции).
SCHEMA_SQL = """
CREATE TABLE general (
    "ID" text PRIMARY KEY,
    number_of_records integer,
    "number_of_SNPs" integer
);

CREATE TABLE snp_sites (
    site_id bigserial PRIMARY KEY,
    site_key text NOT NULL,
    chrom text NOT NULL,
    pos integer NOT NULL,
    ref text NOT NULL,
    alt text NOT NULL,
    qual double precision,
    effect text, impact text, gene text, gene_id text,
    product_name text, strand text, feature_id text, biotype text, rank text,
    hgvs_c text, hgvs_p text, cds_pos text, aa_pos text,
    created_at timestamptz DEFAULT now(),
    CONSTRAINT uq_snp_sites_site_key UNIQUE (site_key)
);

CREATE TABLE sample_snp_alleles (
    sample_id text NOT NULL REFERENCES general("ID") ON DELETE CASCADE,
    site_id bigint NOT NULL REFERENCES snp_sites(site_id) ON DELETE CASCADE,
    allele text NOT NULL,
    dp integer, af double precision, qual double precision, filter text,
    created_at timestamptz DEFAULT now(),
    PRIMARY KEY (sample_id, site_id)
);
CREATE INDEX idx_sample_snp_alleles_sample ON sample_snp_alleles (sample_id);

CREATE TABLE sample_snp_profiles (
    mask_key text NOT NULL,
    sample_id text NOT NULL,
    site_ids bigint[] NOT NULL,
    snp_count integer NOT NULL,
    updated_at timestamptz NOT NULL DEFAULT now(),
    PRIMARY KEY (mask_key, sample_id)
);

CREATE TABLE vcf_table (
    id bigserial PRIMARY KEY,
    sample_id text,
    pos integer,
    alt text
);
CREATE INDEX ix_vcf_table_sample_id ON vcf_table (sample_id);

CREATE TABLE snp_masks (
    mask_key text PRIMARY KEY,
    name text NOT NULL,
    sort_order integer NOT NULL,
    region_count integer NOT NULL DEFAULT 0
);

CREATE TABLE snp_mask_regions (
    mask_key text NOT NULL,
    chrom text NOT NULL,
    pos_range int4range NOT NULL
);
CREATE INDEX idx_snp_mask_regions_range ON snp_mask_regions USING gist (pos_range);
"""

SEED_SQL = """
INSERT INTO general ("ID") VALUES ('ERR4797591'), ('SAMPLE2');

INSERT INTO snp_masks (mask_key, name, sort_order) VALUES
  ('none', 'Без маски', 0),
  ('rlc_lowmap', 'RLC + low-map', 1),
  ('dr_genes', 'DR-гены', 2),
  ('union', 'Объединение', 3);

-- Регион, накрывающий pos 4013: профиль union должен стать короче none.
INSERT INTO snp_mask_regions (mask_key, chrom, pos_range)
VALUES ('union', 'NC_000962.3', int4range(4000, 4100));
"""


def _opted_in():
    """Тест выключен по умолчанию и не запускается вместе с остальными.

    Он поднимает собственный одноразовый Postgres и к рабочим базам не
    подключается, но включать его нужно осознанно:
        TB_LITE_SNP_DB_TEST=1 python3 -m unittest tests.test_import_snp_sql
    """
    if os.environ.get("TB_LITE_SNP_DB_TEST") != "1":
        return False
    if shutil.which("docker") is None:
        return False
    return subprocess.run(["docker", "info"], capture_output=True).returncode == 0


@unittest.skipUnless(_opted_in(), "выключен; включается TB_LITE_SNP_DB_TEST=1")
class ImportSnpSqlIntegrationTest(unittest.TestCase):
    """Прогон import_snp.sql против одноразового Postgres.

    Никогда не подключается к существующим базам: создаёт свой контейнер
    postgres:17 со случайным именем и удаляет его в teardown.
    """

    container: str = ""

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        tmp = Path(cls._tmp.name)

        # Данные готовим тем же путём, что и пайплайн: конвертер -> склейка.
        fixture = REPO_ROOT / "tests" / "data" / "ERR4797591.mini.ann.vcf"
        second = tmp / "SAMPLE2.annotated.ann.vcf"
        lines = [
            line.replace("ERR4797591", "SAMPLE2") if line.startswith("#CHROM") else line
            for line in fixture.read_text(encoding="utf-8").splitlines()
        ]
        second.write_text("\n".join(lines) + "\n", encoding="utf-8")

        shard_dir = tmp / "shards"
        for tag, vcf in (("batch_1", fixture), ("batch_2", second)):
            listing = tmp / f"{tag}.txt"
            listing.write_text(f"{vcf}\n", encoding="utf-8")
            subprocess.run(
                [sys.executable, str(REPO_ROOT / "bin" / "vcf_to_snp_shards.py"),
                 "--file-list", str(listing), "--tag", tag, "-o", str(shard_dir), "-j", "1"],
                check=True,
            )

        cls.data_dir = tmp / "data"
        subprocess.run(
            ["bash", str(REPO_ROOT / "bin" / "merge_snp_shards.sh"),
             "--shard-dir", str(shard_dir), "--out-dir", str(cls.data_dir), "--threads", "1"],
            check=True,
        )

        cls.container = f"tb-lite-snp-test-{uuid.uuid4().hex[:8]}"
        subprocess.run(
            ["docker", "run", "-d", "--rm", "--name", cls.container,
             "-e", "POSTGRES_PASSWORD=test", "-e", "POSTGRES_USER=test",
             "-e", "POSTGRES_DB=test", PG_IMAGE],
            check=True, capture_output=True,
        )
        cls._wait_ready()

        subprocess.run(
            ["docker", "cp", str(cls.data_dir), f"{cls.container}:/data_snp"],
            check=True, capture_output=True,
        )
        cls._psql_script(SCHEMA_SQL)
        cls._psql_script(SEED_SQL)

        cls.sql_path = tmp / "import_snp.sql"
        subprocess.run(
            [sys.executable, str(GENERATOR), "--data-dir", "/data_snp", "-o", str(cls.sql_path)],
            check=True,
        )

    @classmethod
    def tearDownClass(cls):
        if cls.container:
            subprocess.run(["docker", "rm", "-f", cls.container], capture_output=True)
        cls._tmp.cleanup()

    @classmethod
    def _wait_ready(cls, timeout=60):
        deadline = time.time() + timeout
        while time.time() < deadline:
            probe = subprocess.run(
                ["docker", "exec", cls.container, "pg_isready", "-U", "test", "-d", "test"],
                capture_output=True,
            )
            if probe.returncode == 0:
                return
            time.sleep(1)
        raise AssertionError("одноразовый Postgres не поднялся")

    @classmethod
    def _psql(cls, *args, stdin=None):
        return subprocess.run(
            ["docker", "exec", "-i", "-e", "PGPASSWORD=test", cls.container,
             "psql", "-U", "test", "-d", "test", "-v", "ON_ERROR_STOP=1", *args],
            stdin=stdin, capture_output=True, text=stdin is None,
        )

    @classmethod
    def _psql_script(cls, script):
        done = subprocess.run(
            ["docker", "exec", "-i", "-e", "PGPASSWORD=test", cls.container,
             "psql", "-U", "test", "-d", "test", "-v", "ON_ERROR_STOP=1", "-f", "-"],
            input=script, capture_output=True, text=True,
        )
        if done.returncode != 0:
            raise AssertionError(done.stderr)
        return done.stdout

    @classmethod
    def _scalar(cls, statement):
        done = cls._psql("-tAc", statement)
        if done.returncode != 0:
            raise AssertionError(done.stderr)
        return done.stdout.strip()

    def _run_import(self):
        with open(self.sql_path, "rb") as handle:
            done = self._psql("-f", "-", stdin=handle)
        self.assertEqual(done.returncode, 0, done.stderr.decode("utf-8", "replace"))

    def _counts(self):
        return {
            "snp_sites": int(self._scalar("SELECT count(*) FROM snp_sites;")),
            "sample_snp_alleles": int(self._scalar("SELECT count(*) FROM sample_snp_alleles;")),
            "vcf_table": int(self._scalar("SELECT count(*) FROM vcf_table;")),
            "sample_snp_profiles": int(self._scalar("SELECT count(*) FROM sample_snp_profiles;")),
        }

    def test_import_then_reimport_is_idempotent(self):
        self._run_import()
        first = self._counts()

        manifest = dict(
            line.split("\t")
            for line in (self.data_dir / "snp_db_manifest.tsv").read_text().splitlines()
        )
        self.assertEqual(first["snp_sites"], int(manifest["sites_rows"]))
        self.assertEqual(first["sample_snp_alleles"], int(manifest["alleles_rows"]))
        self.assertEqual(first["vcf_table"], int(manifest["vcf_rows"]))
        # 2 образца x 4 маски
        self.assertEqual(first["sample_snp_profiles"], 8)

        self._run_import()
        self.assertEqual(self._counts(), first)

    def test_mask_shrinks_profile(self):
        self._run_import()
        none_count = int(self._scalar(
            "SELECT snp_count FROM sample_snp_profiles "
            "WHERE mask_key='none' AND sample_id='ERR4797591';"
        ))
        union_count = int(self._scalar(
            "SELECT snp_count FROM sample_snp_profiles "
            "WHERE mask_key='union' AND sample_id='ERR4797591';"
        ))
        self.assertGreater(none_count, union_count)

    def test_samples_missing_from_general_are_skipped(self):
        self._run_import()
        self._psql_script("DELETE FROM general WHERE \"ID\" = 'SAMPLE2';")
        self._run_import()
        self.assertEqual(
            self._scalar("SELECT count(*) FROM vcf_table WHERE sample_id='SAMPLE2';"), "0"
        )


if __name__ == "__main__":
    unittest.main()
