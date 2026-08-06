import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "bin" / "filter_table_by_samples.py"


def _run(args):
    return subprocess.run([sys.executable, str(SCRIPT), *args],
                          capture_output=True, text=True)


class FilterTableBySamplesTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)

        self.table = self.tmp / "general.tsv"
        self.table.write_text(
            "ID\tnumber_of_records\tMEAN_COVERAGE\n"
            "ERR1\t1589\t60,76\n"
            "SRR_BAD\t0\t0.0\n"
            "ERR2\t1421\t55,10\n",
            encoding="utf-8",
        )

        self.samples = self.tmp / "passed.txt"
        self.samples.write_text("ERR1\nERR2\n", encoding="utf-8")

        self.out = self.tmp / "out.tsv"

    def tearDown(self):
        self._tmp.cleanup()

    def _rows(self):
        with self.out.open(newline="") as handle:
            return [row for row in csv.reader(handle, delimiter="\t") if row]

    def test_keeps_only_listed_samples(self):
        done = _run(["-i", str(self.table), "-s", str(self.samples),
                     "--id-column", "ID", "-o", str(self.out)])
        self.assertEqual(done.returncode, 0, done.stderr)
        rows = self._rows()
        self.assertEqual(rows[0], ["ID", "number_of_records", "MEAN_COVERAGE"])
        self.assertEqual([row[0] for row in rows[1:]], ["ERR1", "ERR2"])

    def test_preserves_values_verbatim(self):
        _run(["-i", str(self.table), "-s", str(self.samples),
              "--id-column", "ID", "-o", str(self.out)])
        # Десятичная запятая из general.tsv не должна портиться
        self.assertEqual(self._rows()[1][2], "60,76")

    def test_reports_counts_to_stderr(self):
        done = _run(["-i", str(self.table), "-s", str(self.samples),
                     "--id-column", "ID", "-o", str(self.out)])
        self.assertIn("оставлено 2", done.stderr)
        self.assertIn("отброшено 1", done.stderr)

    def test_other_id_column_name(self):
        table = self.tmp / "tbmix.tsv"
        table.write_text("Sample\tLvl1\nERR1\tL2\nSRR_BAD\tL4\n", encoding="utf-8")
        done = _run(["-i", str(table), "-s", str(self.samples),
                     "--id-column", "Sample", "-o", str(self.out)])
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertEqual([row[0] for row in self._rows()[1:]], ["ERR1"])

    def test_missing_column_is_an_error(self):
        done = _run(["-i", str(self.table), "-s", str(self.samples),
                     "--id-column", "sample_id", "-o", str(self.out)])
        self.assertEqual(done.returncode, 1)
        self.assertIn("нет колонки", done.stderr)

    def test_empty_sample_list_keeps_only_header(self):
        empty = self.tmp / "none.txt"
        empty.write_text("", encoding="utf-8")
        done = _run(["-i", str(self.table), "-s", str(empty),
                     "--id-column", "ID", "-o", str(self.out)])
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertEqual(len(self._rows()), 1)


if __name__ == "__main__":
    unittest.main()
