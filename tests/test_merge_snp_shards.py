import csv
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "bin" / "merge_snp_shards.sh"
CONVERTER = REPO_ROOT / "bin" / "vcf_to_snp_shards.py"
FIXTURE = REPO_ROOT / "tests" / "data" / "ERR4797591.mini.ann.vcf"


def _read_gz_rows(path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return [row for row in csv.reader(handle, delimiter="\t") if row]


class MergeSnpShardsTest(unittest.TestCase):
    """Два шарда с пересекающимися сайтами и разным QUAL."""

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        tmp = Path(cls._tmp.name)
        shard_dir = tmp / "shards"

        second = tmp / "SAMPLE2.annotated.ann.vcf"
        lines = []
        for line in FIXTURE.read_text(encoding="utf-8").splitlines():
            if line.startswith("#CHROM"):
                lines.append(line.replace("ERR4797591", "SAMPLE2"))
            elif line.startswith("#"):
                lines.append(line)
            else:
                parts = line.split("\t")
                parts[5] = f"{float(parts[5]) * 3:.2f}"
                lines.append("\t".join(parts))
        second.write_text("\n".join(lines) + "\n", encoding="utf-8")

        for tag, vcf in (("batch_1", FIXTURE), ("batch_2", second)):
            listing = tmp / f"{tag}.txt"
            listing.write_text(f"{vcf}\n", encoding="utf-8")
            subprocess.run(
                [sys.executable, str(CONVERTER), "--file-list", str(listing),
                 "--tag", tag, "-o", str(shard_dir), "-j", "1"],
                check=True,
            )

        cls.outdir = tmp / "out"
        subprocess.run(
            ["bash", str(SCRIPT), "--shard-dir", str(shard_dir),
             "--out-dir", str(cls.outdir), "--threads", "2"],
            check=True,
        )

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_outputs_exist(self):
        for name in ("snp_sites.tsv.gz", "sample_snp_alleles.tsv.gz",
                     "vcf_table.tsv.gz", "samples.txt", "snp_db_manifest.tsv"):
            self.assertTrue((self.outdir / name).exists(), name)

    def test_headers_present(self):
        self.assertEqual(_read_gz_rows(self.outdir / "snp_sites.tsv.gz")[0][0], "site_key")
        self.assertEqual(_read_gz_rows(self.outdir / "sample_snp_alleles.tsv.gz")[0][0], "sample_id")
        self.assertEqual(_read_gz_rows(self.outdir / "vcf_table.tsv.gz")[0], ["sample_id", "pos", "alt"])

    def test_sites_deduplicated_across_shards(self):
        rows = _read_gz_rows(self.outdir / "snp_sites.tsv.gz")[1:]
        self.assertEqual(len(rows), 8)
        keys = [row[0] for row in rows]
        self.assertEqual(len(keys), len(set(keys)))

    def test_sites_keep_max_qual_across_shards(self):
        rows = _read_gz_rows(self.outdir / "snp_sites.tsv.gz")[1:]
        by_pos = {int(row[2]): float(row[5]) for row in rows}
        # У SAMPLE2 QUAL втрое выше, значит должен победить он
        self.assertAlmostEqual(by_pos[1977], 1643.47 * 3, places=1)

    def test_alleles_and_vcf_rows_concatenated(self):
        self.assertEqual(len(_read_gz_rows(self.outdir / "sample_snp_alleles.tsv.gz")) - 1, 16)
        self.assertEqual(len(_read_gz_rows(self.outdir / "vcf_table.tsv.gz")) - 1, 14)

    def test_manifest(self):
        manifest = dict(
            line.split("\t")
            for line in (self.outdir / "snp_db_manifest.tsv").read_text().splitlines()
        )
        self.assertEqual(manifest["shards"], "2")
        self.assertEqual(manifest["samples"], "2")
        self.assertEqual(manifest["sites_rows_before_dedup"], "16")
        self.assertEqual(manifest["sites_rows"], "8")
        self.assertEqual(manifest["alleles_rows"], "16")
        self.assertEqual(manifest["vcf_rows"], "14")


if __name__ == "__main__":
    unittest.main()
