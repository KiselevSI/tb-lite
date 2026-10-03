#!/usr/bin/env python3
"""Сводка производительности batch-прогона TB-Lite (run_batches.sh --benchmark).

Читает трассировку Nextflow (conf/benchmark.config, raw-формат) и тайминги батчей
из <bench-dir> и пишет туда же:

  per_process.tsv  — время/CPU/память по процессам пайплайна
  per_sample.tsv   — ресурсы на образец (+ объём входных данных из fastp и покрытие)
  per_batch.tsv    — wall-clock и пропускная способность по батчам
  summary.tsv      — итоговые метрики (metric, value, unit, description)
  summary.md       — то же в виде текста для Supplementary Table

Только стандартная библиотека: скрипт запускается на хосте, а не в контейнере.
"""

import argparse
import csv
import json
import math
import re
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

RUN_DIR_RE = re.compile(r"^(batch_\d+|final_reports)$")

# Процессы, по которым определяется множество образцов: через них проходит каждый вход.
ANCHOR_PROCESSES = {"FASTP", "SRATOOLS_PREFETCH", "SRA_DETECT_LAYOUT"}
# Скачивание из SRA зависит от сети, поэтому считается отдельно от вычислений.
DOWNLOAD_PROCESSES = {"SRATOOLS_PREFETCH", "SRATOOLS_FASTERQDUMP", "SRA_DETECT_LAYOUT"}
# Образец дошёл до вызова вариантов, т.е. прошёл фильтры покрытия/смеси.
PASSED_QC_PROCESS = "FREEBAYES"

BATCH_LEVEL = "__batch_level__"
MS_PER_HOUR = 3_600_000
GB = 1024 ** 3


# ---------------------------------------------------------------- helpers

def num(value):
    """Число из поля raw-trace; '-', '' и мусор -> None."""
    if value is None:
        return None
    value = value.strip().rstrip("%")
    if value in ("", "-"):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def percentile(values, q):
    """Перцентиль с линейной интерполяцией (как numpy по умолчанию)."""
    if not values:
        return None
    data = sorted(values)
    pos = (len(data) - 1) * q
    lo = math.floor(pos)
    hi = math.ceil(pos)
    return data[lo] + (data[hi] - data[lo]) * (pos - lo)


def median(values):
    return statistics.median(values) if values else None


def fmt(value, digits=3):
    if value is None:
        return "NA"
    if isinstance(value, float):
        if value.is_integer() and abs(value) < 1e15:
            return str(int(value))
        return f"{value:.{digits}f}"
    return str(value)


def short_process(name):
    return name.rsplit(":", 1)[-1]


def normalize_tag(tag):
    """'RD: ERR123' -> 'ERR123'; 'ERR123' -> 'ERR123'."""
    tag = (tag or "").strip()
    if ": " in tag:
        tag = tag.split(": ", 1)[1].strip()
    return tag


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, columns, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            writer.writerow([fmt(row.get(col)) for col in columns])


# ---------------------------------------------------------------- loading

class Task:
    __slots__ = (
        "run", "process", "sample_tag", "status", "hash", "cpus", "realtime",
        "pcpu", "peak_rss", "start", "complete", "read_bytes", "write_bytes",
        "rchar", "wchar", "from_attempt",
    )

    def __init__(self, run, row):
        self.run = run
        self.process = short_process(row.get("process", ""))
        self.sample_tag = normalize_tag(row.get("tag"))
        self.status = (row.get("status") or "").strip()
        self.hash = (row.get("hash") or "").strip()
        self.from_attempt = False
        self._load_metrics(row)

    def _load_metrics(self, row):
        self.cpus = num(row.get("cpus")) or 1.0
        self.realtime = num(row.get("realtime"))
        self.pcpu = num(row.get("%cpu"))
        self.peak_rss = num(row.get("peak_rss"))
        self.start = num(row.get("start"))
        self.complete = num(row.get("complete"))
        self.read_bytes = num(row.get("read_bytes"))
        self.write_bytes = num(row.get("write_bytes"))
        self.rchar = num(row.get("rchar"))
        self.wchar = num(row.get("wchar"))

    def take_metrics_from(self, row):
        self._load_metrics(row)
        self.from_attempt = True

    @property
    def has_metrics(self):
        return self.realtime is not None

    @property
    def cpu_hours(self):
        """Фактически использованное CPU-время (realtime × %cpu)."""
        if self.realtime is None or self.pcpu is None:
            return None
        return self.realtime * self.pcpu / 100.0 / MS_PER_HOUR

    @property
    def alloc_core_hours(self):
        """Выделенные ядро-часы (realtime × cpus) — то, что «занято» в планировщике."""
        if self.realtime is None:
            return None
        return self.realtime * self.cpus / MS_PER_HOUR


def load_traces(bench_dir):
    """Основные trace (последняя попытка каждого запуска) + индекс прошлых попыток."""
    tasks = []
    for run_dir in sorted(p for p in bench_dir.iterdir() if p.is_dir() and RUN_DIR_RE.match(p.name)):
        trace = run_dir / "trace.tsv"
        if not trace.is_file():
            continue
        tasks.extend(Task(run_dir.name, row) for row in read_tsv(trace))

    attempt_completed = {}
    attempt_failed = []
    attempts_dir = bench_dir / "attempts"
    if attempts_dir.is_dir():
        for trace in sorted(attempts_dir.glob("*/trace.tsv")):
            run = trace.parent.name.split(".", 1)[0]
            for row in read_tsv(trace):
                status = (row.get("status") or "").strip()
                if status == "COMPLETED" and row.get("hash"):
                    attempt_completed[row["hash"].strip()] = row
                elif status in ("FAILED", "ABORTED"):
                    attempt_failed.append(Task(run, row))
    return tasks, attempt_completed, attempt_failed


def load_runs(bench_dir):
    path = bench_dir / "batches.tsv"
    if not path.is_file():
        return []
    runs = []
    for row in read_tsv(path):
        runs.append({
            "run": row["run"],
            "samples": num(row.get("samples")),
            "wall_sec": num(row.get("wall_sec")),
            "exit_code": num(row.get("exit_code")),
            "work_bytes": num(row.get("work_bytes")),
        })
    return runs


def load_environment(bench_dir):
    path = bench_dir / "environment.tsv"
    if not path.is_file():
        return {}
    return {row["key"]: row["value"] for row in read_tsv(path)}


def load_fastp(outdir):
    data = {}
    if outdir is None:
        return data
    for path in outdir.glob("fastp/*/*.fastp.json"):
        try:
            before = json.loads(path.read_text())["summary"]["before_filtering"]
        except (OSError, ValueError, KeyError):
            continue
        data[path.parent.name] = {
            "reads": float(before.get("total_reads", 0)),
            "bases": float(before.get("total_bases", 0)),
        }
    return data


def load_coverage(outdir):
    data = {}
    if outdir is None:
        return data
    for path in outdir.glob("stats/picard/wgs/*/*coverage_metrics"):
        try:
            lines = path.read_text().splitlines()
        except OSError:
            continue
        for i, line in enumerate(lines[:-1]):
            if line.startswith("GENOME_TERRITORY"):
                header = line.split("\t")
                values = lines[i + 1].split("\t")
                if "MEAN_COVERAGE" in header:
                    data[path.parent.name] = num(values[header.index("MEAN_COVERAGE")])
                break
    return data


def dir_size(path):
    if path is None or not path.is_dir():
        return None
    try:
        out = subprocess.run(["du", "-sb", str(path)], capture_output=True, text=True, check=True)
        return float(out.stdout.split()[0])
    except (OSError, subprocess.CalledProcessError, ValueError, IndexError):
        return None


# ---------------------------------------------------------------- aggregation

def resolve_cached(tasks, attempt_completed):
    """CACHED-задачи берут метрики из прошлой попытки с тем же hash, если она есть."""
    missing = 0
    for task in tasks:
        if task.status != "CACHED":
            continue
        row = attempt_completed.get(task.hash)
        if row is not None:
            task.take_metrics_from(row)
        elif not task.has_metrics:
            missing += 1
    return missing


def is_measured(task):
    return task.status in ("COMPLETED", "CACHED") and task.has_metrics


def sum_or_none(values):
    values = [v for v in values if v is not None]
    return sum(values) if values else None


def aggregate_processes(tasks, failed_tasks, total_cpu_h):
    by_process = defaultdict(list)
    for task in tasks:
        if is_measured(task):
            by_process[task.process].append(task)
    failed = defaultdict(int)
    for task in failed_tasks:
        failed[task.process] += 1

    rows = []
    for process in sorted(set(by_process) | set(failed)):
        items = by_process.get(process, [])
        realtime_min = [t.realtime / 60000 for t in items]
        rss_gb = [t.peak_rss / GB for t in items if t.peak_rss is not None]
        pcpu = [t.pcpu for t in items if t.pcpu is not None]
        cpu_h = sum_or_none(t.cpu_hours for t in items) or 0.0
        rows.append({
            "process": process,
            "n_tasks": len(items),
            "n_failed": failed.get(process, 0),
            "cpus_allocated": max((t.cpus for t in items), default=None),
            "realtime_median_min": median(realtime_min),
            "realtime_p95_min": percentile(realtime_min, 0.95),
            "realtime_max_min": max(realtime_min, default=None),
            "peak_rss_median_gb": median(rss_gb),
            "peak_rss_max_gb": max(rss_gb, default=None),
            "pcpu_mean": statistics.fmean(pcpu) if pcpu else None,
            "cpu_hours": cpu_h,
            "alloc_core_hours": sum_or_none(t.alloc_core_hours for t in items),
            "cpu_share_pct": 100.0 * cpu_h / total_cpu_h if total_cpu_h else None,
        })
    rows.sort(key=lambda r: r["cpu_hours"] or 0.0, reverse=True)
    return rows


def aggregate_samples(tasks, samples, fastp, coverage):
    by_sample = defaultdict(list)
    for task in tasks:
        if task.sample_tag in samples and is_measured(task):
            by_sample[task.sample_tag].append(task)

    rows = []
    for sample in sorted(samples):
        items = by_sample.get(sample, [])
        compute = [t for t in items if t.process not in DOWNLOAD_PROCESSES]
        download = [t for t in items if t.process in DOWNLOAD_PROCESSES]
        starts = [t.start for t in compute if t.start is not None]
        ends = [t.complete for t in compute if t.complete is not None]
        all_starts = [t.start for t in items if t.start is not None]
        all_ends = [t.complete for t in items if t.complete is not None]
        rss = [(t.peak_rss, t.process) for t in items if t.peak_rss is not None]
        top_rss = max(rss, default=(None, None))
        inp = fastp.get(sample, {})
        rows.append({
            "sample": sample,
            "batch": samples[sample],
            "passed_qc": int(any(t.process == PASSED_QC_PROCESS for t in items)),
            "n_tasks": len(items),
            "cpu_hours": sum_or_none(t.cpu_hours for t in compute),
            "alloc_core_hours": sum_or_none(t.alloc_core_hours for t in compute),
            "sum_realtime_h": sum_or_none(t.realtime / MS_PER_HOUR for t in compute),
            "latency_h": (max(ends) - min(starts)) / MS_PER_HOUR if starts and ends else None,
            "latency_with_download_h": (max(all_ends) - min(all_starts)) / MS_PER_HOUR if all_starts and all_ends else None,
            "download_h": sum_or_none(t.realtime / MS_PER_HOUR for t in download),
            "peak_rss_gb": top_rss[0] / GB if top_rss[0] is not None else None,
            "peak_rss_process": top_rss[1],
            "read_gb": sum_or_none((t.rchar or 0) / GB for t in compute),
            "write_gb": sum_or_none((t.wchar or 0) / GB for t in compute),
            "input_reads": inp.get("reads"),
            "input_gbases": inp["bases"] / 1e9 if "bases" in inp else None,
            "mean_coverage": coverage.get(sample),
        })
    return rows


def aggregate_batches(runs, tasks):
    cpu_by_run = defaultdict(float)
    failed_by_run = defaultdict(int)
    for task in tasks:
        if is_measured(task) and task.cpu_hours is not None:
            cpu_by_run[task.run] += task.cpu_hours
        if task.status in ("FAILED", "ABORTED"):
            failed_by_run[task.run] += 1

    rows = []
    for run in runs:
        wall_h = run["wall_sec"] / 3600 if run["wall_sec"] else None
        ok = run["exit_code"] == 0
        samples_per_h = run["samples"] / wall_h if ok and wall_h and run["run"] != "final_reports" else None
        rows.append({
            "run": run["run"],
            "exit_code": run["exit_code"],
            "samples": run["samples"],
            "wall_h": wall_h,
            "samples_per_hour": samples_per_h,
            # CPU-часы есть только у последней попытки запуска (trace перезаписывается)
            "cpu_hours": cpu_by_run.get(run["run"]) if ok else None,
            "work_gb": run["work_bytes"] / GB if run["work_bytes"] is not None else None,
            "failed_tasks": failed_by_run.get(run["run"], 0) if ok else None,
        })
    return rows


def find_samples(tasks):
    """Образец -> батч. Образцы = теги задач «якорных» процессов (вход каждого образца)."""
    samples = {}
    for task in tasks:
        if task.process in ANCHOR_PROCESSES and task.sample_tag:
            samples.setdefault(task.sample_tag, task.run)
    return samples


# ---------------------------------------------------------------- summary

def build_summary(env, runs, tasks, failed_tasks, samples, sample_rows, batch_rows,
                  cached_missing, outdir_bytes):
    metrics = []

    def add(key, value, unit, description):
        metrics.append({"metric": key, "value": value, "unit": unit, "description": description})

    n_samples = len(samples)
    n_passed = sum(r["passed_qc"] for r in sample_rows)
    wall_total_h = sum((r["wall_sec"] or 0) for r in runs) / 3600 if runs else None
    batch_speeds = [r["samples_per_hour"] for r in batch_rows if r["samples_per_hour"] is not None]

    measured = [t for t in tasks if is_measured(t)]
    compute = [t for t in measured if t.process not in DOWNLOAD_PROCESSES]
    total_cpu_h = sum_or_none(t.cpu_hours for t in measured)
    total_alloc_h = sum_or_none(t.alloc_core_hours for t in measured)
    batch_level_cpu_h = sum_or_none(t.cpu_hours for t in compute if t.sample_tag not in samples)

    per_sample_cpu = [r["cpu_hours"] for r in sample_rows if r["cpu_hours"] is not None]
    per_sample_latency = [r["latency_h"] for r in sample_rows if r["latency_h"] is not None]
    per_sample_download = [r["download_h"] for r in sample_rows if r["download_h"] is not None]
    per_sample_rss = [r["peak_rss_gb"] for r in sample_rows if r["peak_rss_gb"] is not None]
    gbases = [r["input_gbases"] for r in sample_rows if r["input_gbases"] is not None]
    reads = [r["input_reads"] for r in sample_rows if r["input_reads"] is not None]
    cov = [r["mean_coverage"] for r in sample_rows if r["mean_coverage"] is not None]
    work_gb = [r["work_gb"] for r in batch_rows if r["work_gb"] is not None and r["run"] != "final_reports"]

    peak = max((t for t in measured if t.peak_rss is not None), key=lambda t: t.peak_rss, default=None)
    executor_cpus = num(env.get("executor_cpus"))

    n_failed = sum(1 for t in tasks if t.status in ("FAILED", "ABORTED")) + len(failed_tasks)
    n_cached = sum(1 for t in tasks if t.status == "CACHED")
    n_total = len(tasks) + len(failed_tasks)

    add("samples_input", n_samples, "samples", "Образцов на входе (прошли FASTP/загрузку SRA)")
    add("samples_passed_qc", n_passed, "samples", f"Образцов, дошедших до {PASSED_QC_PROCESS} (прошли фильтры)")
    add("batches", sum(1 for r in batch_rows if r["run"] != "final_reports" and r["exit_code"] == 0),
        "batches", "Успешно завершённых батчей")
    add("wall_clock_total_h", wall_total_h, "h",
        "Суммарное wall-clock время всех запусков, включая упавшие попытки и final_reports")
    add("throughput_samples_per_hour", n_samples / wall_total_h if wall_total_h else None, "samples/h",
        "End-to-end пропускная способность: samples_input / wall_clock_total_h")
    add("throughput_samples_per_day", 24 * n_samples / wall_total_h if wall_total_h else None, "samples/day",
        "То же в пересчёте на сутки")
    add("throughput_batch_median_samples_per_hour", median(batch_speeds), "samples/h",
        "Медиана пропускной способности по успешным батчам (без final_reports)")
    add("gbases_per_hour", sum(gbases) / wall_total_h if gbases and wall_total_h else None, "Gbp/h",
        "Обработано сырых оснований (fastp, до фильтрации) в час wall-clock")
    add("cpu_hours_total", total_cpu_h, "CPU·h", "Использованное CPU-время всех задач (realtime × %cpu)")
    add("alloc_core_hours_total", total_alloc_h, "core·h", "Выделенные ядро-часы (realtime × cpus)")
    add("cpu_efficiency_pct", 100 * total_cpu_h / total_alloc_h if total_cpu_h and total_alloc_h else None, "%",
        "Использованное CPU / выделенное CPU")
    add("executor_cpus", executor_cpus, "cores", "Лимит ядер executor.cpus")
    add("node_utilization_pct",
        100 * total_cpu_h / (wall_total_h * executor_cpus) if total_cpu_h and wall_total_h and executor_cpus else None,
        "%", "CPU-часы / (wall-clock × executor.cpus)")
    add("cpu_hours_per_sample_median", median(per_sample_cpu), "CPU·h",
        "Медиана CPU-часов на образец (без загрузки SRA и общих задач батча)")
    add("cpu_hours_per_sample_p25", percentile(per_sample_cpu, 0.25), "CPU·h", "25-й перцентиль")
    add("cpu_hours_per_sample_p75", percentile(per_sample_cpu, 0.75), "CPU·h", "75-й перцентиль")
    add("cpu_hours_per_sample_mean", statistics.fmean(per_sample_cpu) if per_sample_cpu else None, "CPU·h",
        "Среднее CPU-часов на образец")
    add("batch_level_cpu_hours_per_sample",
        batch_level_cpu_h / n_samples if batch_level_cpu_h is not None and n_samples else None, "CPU·h",
        "Общие задачи батча (индексы, итоговые таблицы) в расчёте на образец")
    add("latency_per_sample_median_h", median(per_sample_latency), "h",
        "Медиана времени от первой до последней задачи образца (без загрузки SRA)")
    add("latency_per_sample_p95_h", percentile(per_sample_latency, 0.95), "h", "95-й перцентиль")
    add("sra_download_per_sample_median_h", median(per_sample_download), "h",
        "Медиана времени загрузки/распаковки SRA на образец (0 образцов — вход FASTQ)")
    add("peak_rss_per_sample_median_gb", median(per_sample_rss), "GB",
        "Медиана максимального peak RSS среди задач образца")
    add("peak_rss_max_gb", peak.peak_rss / GB if peak else None, "GB",
        f"Максимальный peak RSS одной задачи ({peak.process if peak else 'NA'})")
    add("input_reads_per_sample_median", median(reads), "reads", "Медиана числа сырых чтений (fastp)")
    add("input_gbases_per_sample_median", median(gbases), "Gbp", "Медиана объёма сырых данных на образец (fastp)")
    add("mean_coverage_median", median(cov), "x", "Медиана среднего покрытия (Picard CollectWgsMetrics)")
    add("work_dir_peak_gb_max", max(work_gb, default=None), "GB",
        "Максимальный объём work/ за батч (перед очисткой)")
    add("outdir_gb_per_sample", outdir_bytes / GB / n_samples if outdir_bytes and n_samples else None, "GB",
        "Объём outdir в расчёте на образец")
    add("tasks_total", n_total, "tasks", "Всего задач Nextflow (включая упавшие попытки)")
    add("tasks_failed", n_failed, "tasks", "Задач со статусом FAILED/ABORTED (ретраи и падения)")
    add("task_failure_rate_pct", 100 * n_failed / n_total if n_total else None, "%", "Доля упавших задач")
    add("tasks_cached", n_cached, "tasks", "Задач, взятых из кэша (-resume)")

    warnings = []
    if cached_missing:
        warnings.append(
            f"{cached_missing} CACHED-задач без метрик исключены из расчёта: для статьи "
            "используйте чистый прогон без -resume."
        )
    if any(r["exit_code"] not in (0, None) for r in runs):
        warnings.append("Есть упавшие попытки запусков: их wall-clock входит в wall_clock_total_h.")
    if not runs:
        warnings.append("batches.tsv не найден — wall-clock и samples/h не посчитаны.")
    return metrics, warnings


def write_markdown(path, env, metrics, process_rows, warnings):
    lines = ["# TB-Lite benchmark", ""]
    if env:
        lines += ["## Окружение", "", "| параметр | значение |", "|---|---|"]
        lines += [f"| {k} | {v} |" for k, v in env.items()]
        lines.append("")
    lines += ["## Итоговые метрики", "", "| метрика | значение | ед. | описание |", "|---|---|---|---|"]
    lines += [f"| {m['metric']} | {fmt(m['value'])} | {m['unit']} | {m['description']} |" for m in metrics]
    lines.append("")
    lines += [
        "## Процессы (по убыванию CPU-часов)", "",
        "| процесс | задач | cpus | realtime медиана, мин | p95, мин | peak RSS макс, GB | CPU·h | доля CPU, % |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in process_rows:
        lines.append(
            f"| {r['process']} | {r['n_tasks']} | {fmt(r['cpus_allocated'])} | {fmt(r['realtime_median_min'], 2)} "
            f"| {fmt(r['realtime_p95_min'], 2)} | {fmt(r['peak_rss_max_gb'], 2)} | {fmt(r['cpu_hours'], 3)} "
            f"| {fmt(r['cpu_share_pct'], 1)} |"
        )
    if warnings:
        lines += ["", "## Предупреждения", ""]
        lines += [f"- {w}" for w in warnings]
    path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------- main

PROCESS_COLUMNS = [
    "process", "n_tasks", "n_failed", "cpus_allocated", "realtime_median_min", "realtime_p95_min",
    "realtime_max_min", "peak_rss_median_gb", "peak_rss_max_gb", "pcpu_mean", "cpu_hours",
    "alloc_core_hours", "cpu_share_pct",
]
SAMPLE_COLUMNS = [
    "sample", "batch", "passed_qc", "n_tasks", "cpu_hours", "alloc_core_hours", "sum_realtime_h",
    "latency_h", "latency_with_download_h", "download_h", "peak_rss_gb", "peak_rss_process",
    "read_gb", "write_gb", "input_reads", "input_gbases", "mean_coverage",
]
BATCH_COLUMNS = [
    "run", "exit_code", "samples", "wall_h", "samples_per_hour", "cpu_hours", "work_gb", "failed_tasks",
]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bench-dir", required=True, type=Path, help="Каталог <outdir>/benchmark")
    parser.add_argument("--outdir", type=Path, default=None,
                        help="Каталог результатов прогона (для fastp/Picard и объёма outdir)")
    parser.add_argument("--no-du", action="store_true", help="Не считать объём outdir (du на больших прогонах долгий)")
    args = parser.parse_args(argv)

    bench_dir = args.bench_dir
    if not bench_dir.is_dir():
        parser.error(f"нет каталога {bench_dir}")

    tasks, attempt_completed, attempt_failed = load_traces(bench_dir)
    if not tasks:
        print(f"Ошибка: в {bench_dir} нет trace.tsv (batch_*/, final_reports/)", file=sys.stderr)
        return 1

    cached_missing = resolve_cached(tasks, attempt_completed)
    runs = load_runs(bench_dir)
    env = load_environment(bench_dir)
    samples = find_samples(tasks)
    fastp = load_fastp(args.outdir)
    coverage = load_coverage(args.outdir)
    outdir_bytes = None if args.no_du else dir_size(args.outdir)

    total_cpu_h = sum_or_none(t.cpu_hours for t in tasks if is_measured(t))
    failed_all = [t for t in tasks if t.status in ("FAILED", "ABORTED")] + attempt_failed

    process_rows = aggregate_processes(tasks, failed_all, total_cpu_h)
    sample_rows = aggregate_samples(tasks, samples, fastp, coverage)
    batch_rows = aggregate_batches(runs, tasks)
    metrics, warnings = build_summary(
        env, runs, tasks, attempt_failed, samples, sample_rows, batch_rows,
        cached_missing, outdir_bytes,
    )

    write_tsv(bench_dir / "per_process.tsv", PROCESS_COLUMNS, process_rows)
    write_tsv(bench_dir / "per_sample.tsv", SAMPLE_COLUMNS, sample_rows)
    write_tsv(bench_dir / "per_batch.tsv", BATCH_COLUMNS, batch_rows)
    write_tsv(bench_dir / "summary.tsv", ["metric", "value", "unit", "description"], metrics)
    write_markdown(bench_dir / "summary.md", env, metrics, process_rows, warnings)

    for warning in warnings:
        print(f"ПРЕДУПРЕЖДЕНИЕ: {warning}", file=sys.stderr)
    print(f"Образцов: {len(samples)}; сводка: {bench_dir / 'summary.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
