import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "bin" / "benchmark_summary.py"

TRACE_COLUMNS = [
    "task_id", "hash", "native_id", "process", "tag", "name", "status", "exit", "attempt", "cpus",
    "memory", "submit", "start", "complete", "duration", "realtime", "%cpu", "%mem", "peak_rss",
    "peak_vmem", "rchar", "wchar", "read_bytes", "write_bytes", "hostname", "cpu_model",
]
HOUR_MS = 3_600_000
GB = 1024 ** 3


def _task(hash_, process, tag, status, cpus, start_h, realtime_h, pcpu, rss_gb):
    start = int(start_h * HOUR_MS)
    realtime = int(realtime_h * HOUR_MS)
    return {
        "task_id": "1", "hash": hash_, "native_id": "1", "process": f"TBLITE:SUB:{process}",
        "tag": tag, "name": process, "status": status, "exit": "0", "attempt": "1",
        "cpus": str(cpus), "memory": "-", "submit": str(start), "start": str(start),
        "complete": str(start + realtime), "duration": str(realtime), "realtime": str(realtime),
        "%cpu": str(pcpu), "%mem": "1.0", "peak_rss": str(int(rss_gb * GB)), "peak_vmem": "-",
        "rchar": str(GB), "wchar": str(GB), "read_bytes": "0", "write_bytes": "0",
        "hostname": "host", "cpu_model": "cpu",
    }


def _write_tsv(path, columns, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class BenchmarkSummaryTest(unittest.TestCase):
    """Два батча (второй — после падения и -resume) и final_reports."""

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        tmp = Path(cls._tmp.name)
        outdir = tmp / "results"
        bench = outdir / "benchmark"

        # batch_1: S1 проходит фильтры, S2 отсеян после FASTP; плюс общая задача батча.
        _write_tsv(bench / "batch_1" / "trace.tsv", TRACE_COLUMNS, [
            _task("aa/01", "FASTP", "S1", "COMPLETED", 4, 0.0, 0.5, 200, 1.0),
            _task("aa/02", "RD", "RD: S1", "COMPLETED", 1, 0.5, 0.5, 100, 3.0),
            _task("aa/03", "FREEBAYES", "S1", "COMPLETED", 1, 0.5, 1.0, 100, 2.0),
            _task("aa/04", "FASTP", "S2", "COMPLETED", 4, 0.0, 0.25, 400, 0.5),
            _task("aa/05", "FINAL_TABLE", "Final Table", "COMPLETED", 1, 1.5, 1.0, 100, 0.5),
        ])
        # batch_2, первая попытка: S3 FASTP посчитан, FREEBAYES упал.
        _write_tsv(bench / "attempts" / "batch_2.1700000000" / "trace.tsv", TRACE_COLUMNS, [
            _task("bb/01", "FASTP", "S3", "COMPLETED", 4, 0.0, 1.0, 100, 1.0),
            _task("bb/02", "FREEBAYES", "S3", "FAILED", 1, 1.0, 0.1, 100, 1.0),
        ])
        # batch_2 после -resume: FASTP из кэша (метрики берутся из attempts), плюс
        # CACHED-задача без метрик и без пары в attempts.
        cached = _task("bb/01", "FASTP", "S3", "CACHED", 4, 0.0, 0.0, 0, 0.0)
        for field in ("realtime", "%cpu", "peak_rss", "start", "complete"):
            cached[field] = "-"
        orphan = _task("bb/09", "RD", "RD: S3", "CACHED", 1, 0.0, 0.0, 0, 0.0)
        for field in ("realtime", "%cpu", "peak_rss", "start", "complete"):
            orphan[field] = "-"
        _write_tsv(bench / "batch_2" / "trace.tsv", TRACE_COLUMNS, [
            cached,
            orphan,
            _task("bb/03", "FREEBAYES", "S3", "COMPLETED", 1, 2.0, 1.0, 100, 2.0),
        ])
        _write_tsv(bench / "final_reports" / "trace.tsv", TRACE_COLUMNS, [
            _task("cc/01", "FINAL_TABLE", "Final Table", "COMPLETED", 1, 5.0, 0.5, 100, 4.0),
        ])
        _write_tsv(bench / "batches.tsv",
                   ["run", "samples", "start_epoch", "end_epoch", "wall_sec", "exit_code", "work_bytes"], [
            {"run": "batch_1", "samples": "2", "start_epoch": "0", "end_epoch": "7200",
             "wall_sec": "7200", "exit_code": "0", "work_bytes": str(10 * GB)},
            {"run": "batch_2", "samples": "1", "start_epoch": "7200", "end_epoch": "10800",
             "wall_sec": "3600", "exit_code": "1", "work_bytes": str(5 * GB)},
            {"run": "batch_2", "samples": "1", "start_epoch": "10800", "end_epoch": "14400",
             "wall_sec": "3600", "exit_code": "0", "work_bytes": str(6 * GB)},
            {"run": "final_reports", "samples": "3", "start_epoch": "14400", "end_epoch": "18000",
             "wall_sec": "3600", "exit_code": "0", "work_bytes": "NA"},
        ])
        _write_tsv(bench / "environment.tsv", ["key", "value"], [
            {"key": "executor_cpus", "value": "4"},
            {"key": "cpu_model", "value": "Test CPU"},
        ])
        fastp = outdir / "fastp" / "S1" / "S1.fastp.json"
        fastp.parent.mkdir(parents=True)
        fastp.write_text(json.dumps({"summary": {"before_filtering": {
            "total_reads": 1_000_000, "total_bases": 2_000_000_000}}}))

        cls.result = subprocess.run(
            [sys.executable, str(SCRIPT), "--bench-dir", str(bench), "--outdir", str(outdir), "--no-du"],
            capture_output=True, text=True,
        )
        cls.bench = bench

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def setUp(self):
        self.assertEqual(self.result.returncode, 0, self.result.stderr)

    def _summary(self):
        return {r["metric"]: r["value"] for r in _read_tsv(self.bench / "summary.tsv")}

    def test_samples_found_from_anchor_processes(self):
        rows = {r["sample"]: r for r in _read_tsv(self.bench / "per_sample.tsv")}
        self.assertEqual(sorted(rows), ["S1", "S2", "S3"])
        self.assertEqual(rows["S1"]["batch"], "batch_1")
        self.assertEqual(rows["S3"]["batch"], "batch_2")

    def test_prefixed_tag_counts_towards_sample(self):
        s1 = {r["sample"]: r for r in _read_tsv(self.bench / "per_sample.tsv")}["S1"]
        self.assertEqual(s1["n_tasks"], "3")
        # 0.5×2 + 0.5×1 + 1.0×1 CPU·h
        self.assertAlmostEqual(float(s1["cpu_hours"]), 2.5)
        self.assertAlmostEqual(float(s1["latency_h"]), 1.5)
        self.assertAlmostEqual(float(s1["peak_rss_gb"]), 3.0)
        self.assertEqual(s1["peak_rss_process"], "RD")
        self.assertAlmostEqual(float(s1["input_gbases"]), 2.0)

    def test_passed_qc_flag(self):
        rows = {r["sample"]: r["passed_qc"] for r in _read_tsv(self.bench / "per_sample.tsv")}
        self.assertEqual(rows, {"S1": "1", "S2": "0", "S3": "1"})

    def test_cached_task_takes_metrics_from_previous_attempt(self):
        s3 = {r["sample"]: r for r in _read_tsv(self.bench / "per_sample.tsv")}["S3"]
        # FASTP из attempts (1.0 CPU·h) + FREEBAYES после resume (1.0); orphan RD без метрик пропущен
        self.assertAlmostEqual(float(s3["cpu_hours"]), 2.0)
        self.assertEqual(s3["n_tasks"], "2")
        self.assertIn("1 CACHED", (self.bench / "summary.md").read_text())

    def test_batch_level_tasks_not_assigned_to_samples(self):
        summary = self._summary()
        # FINAL_TABLE: 1.0 + 0.5 CPU·h на 3 образца
        self.assertAlmostEqual(float(summary["batch_level_cpu_hours_per_sample"]), 0.5)

    def test_throughput_uses_all_wall_clock(self):
        summary = self._summary()
        self.assertEqual(summary["samples_input"], "3")
        self.assertEqual(summary["samples_passed_qc"], "2")
        self.assertAlmostEqual(float(summary["wall_clock_total_h"]), 5.0)
        self.assertAlmostEqual(float(summary["throughput_samples_per_hour"]), 0.6)
        # батчи: 2 обр / 2 ч = 1.0 и 1 обр / 1 ч = 1.0
        self.assertAlmostEqual(float(summary["throughput_batch_median_samples_per_hour"]), 1.0)
        self.assertEqual(summary["tasks_failed"], "1")
        self.assertAlmostEqual(float(summary["work_dir_peak_gb_max"]), 10.0)

    def test_per_process_sorted_by_cpu(self):
        rows = _read_tsv(self.bench / "per_process.tsv")
        self.assertEqual(rows[0]["process"], "FASTP")
        freebayes = next(r for r in rows if r["process"] == "FREEBAYES")
        self.assertEqual(freebayes["n_failed"], "1")

    def test_failed_batch_attempt_has_no_speed(self):
        rows = _read_tsv(self.bench / "per_batch.tsv")
        failed = [r for r in rows if r["exit_code"] == "1"]
        self.assertEqual(len(failed), 1)
        self.assertEqual(failed[0]["samples_per_hour"], "NA")


if __name__ == "__main__":
    unittest.main()
