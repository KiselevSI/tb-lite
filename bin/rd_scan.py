#!/usr/bin/env python3
"""
Memory-bounded, resumable batch detection of known and novel RD deletions
from mosdepth *.per-base.bed[.gz] files.

Key properties
--------------
* -i accepts either one mosdepth file or a directory.
* Coverage files are streamed; pandas/numpy are not used.
* Each sample is read twice when mosdepth.summary.txt is unavailable:
  pass 1 calculates chromosome-wide weighted mean coverage;
  pass 2 evaluates known RDs and scans for novel low-coverage intervals.
* Parallelism is bounded: at most --max-pending samples are queued.
* Results and progress are committed transactionally to SQLite, allowing a
  stopped run to continue safely without duplicating completed samples.
* The final output is one TSV containing known and novel deletion calls.
* Progress reporting includes percentage, rolling and average speed, elapsed
  time, ETA, queue size, successful samples, errors, and deletion-row count.

Coordinate convention
---------------------
Coordinates are preserved exactly as supplied by RD.bed and mosdepth BED files
and are treated internally as BED-style 0-based, half-open intervals [start,end).
"""

from __future__ import annotations

import argparse
import csv
from collections import deque
import gzip
import hashlib
import json
import logging
import math
import os
import sqlite3
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


LOG = logging.getLogger("rd_scan_batch")


@dataclass(frozen=True)
class KnownRD:
    chrom: str
    start: int
    end: int
    name: str

    @property
    def size(self) -> int:
        return self.end - self.start


@dataclass
class ChromStats:
    bases: int = 0
    weighted_coverage: float = 0.0
    min_start: Optional[int] = None
    max_end: Optional[int] = None
    records: int = 0

    def add(self, start: int, end: int, coverage: float) -> None:
        length = end - start
        self.bases += length
        self.weighted_coverage += length * coverage
        self.records += 1
        if self.min_start is None or start < self.min_start:
            self.min_start = start
        if self.max_end is None or end > self.max_end:
            self.max_end = end

    @property
    def mean(self) -> float:
        return self.weighted_coverage / self.bases if self.bases > 0 else 0.0


@dataclass
class KnownAccumulator:
    covered_bases: int = 0
    weighted_coverage: float = 0.0
    low_bases: int = 0
    max_low_run: int = 0
    current_low_run: int = 0
    last_low_end: Optional[int] = None

    def add(self, start: int, end: int, coverage: float, low_threshold: float) -> None:
        length = end - start
        if length <= 0:
            return
        self.covered_bases += length
        self.weighted_coverage += length * coverage

        if coverage < low_threshold:
            self.low_bases += length
            if self.last_low_end == start:
                self.current_low_run += length
            else:
                self.current_low_run = length
            self.last_low_end = end
            if self.current_low_run > self.max_low_run:
                self.max_low_run = self.current_low_run
        else:
            self.current_low_run = 0
            self.last_low_end = None


@dataclass
class Candidate:
    chrom: str
    start: int
    last_low_end: int
    weighted_coverage: float
    low_bases: int
    max_low_run: int
    current_low_run: int
    pending_gap_bases: int = 0
    pending_gap_weighted: float = 0.0

    @property
    def end(self) -> int:
        return self.last_low_end

    @property
    def size(self) -> int:
        return self.last_low_end - self.start

    def add_low(self, start: int, end: int, coverage: float) -> None:
        length = end - start
        if self.pending_gap_bases:
            self.weighted_coverage += self.pending_gap_weighted
            self.current_low_run = 0
            self.pending_gap_bases = 0
            self.pending_gap_weighted = 0.0

        self.weighted_coverage += length * coverage
        self.low_bases += length
        if self.last_low_end == start and self.current_low_run > 0:
            self.current_low_run += length
        else:
            self.current_low_run = length
        self.last_low_end = end
        if self.current_low_run > self.max_low_run:
            self.max_low_run = self.current_low_run

    def add_gap(self, length: int, coverage: float) -> None:
        self.pending_gap_bases += length
        self.pending_gap_weighted += length * coverage


@dataclass
class NovelDetector:
    chrom: str
    baseline: float
    low_factor: float
    merge_distance: int
    min_size: int
    max_size: int
    candidate: Optional[Candidate] = None
    last_end: Optional[int] = None

    @property
    def low_threshold(self) -> float:
        return self.baseline * self.low_factor

    def feed(self, start: int, end: int, coverage: float, output: List[dict]) -> None:
        if end <= start:
            return

        is_low = coverage < self.low_threshold
        length = end - start

        if self.candidate is None:
            if is_low:
                self.candidate = Candidate(
                    chrom=self.chrom,
                    start=start,
                    last_low_end=end,
                    weighted_coverage=length * coverage,
                    low_bases=length,
                    max_low_run=length,
                    current_low_run=length,
                )
            return

        if is_low:
            # A candidate is finalized as soon as a normal-coverage gap exceeds
            # merge_distance, so any remaining pending gap is safe to merge.
            self.candidate.add_low(start, end, coverage)
            return

        self.candidate.add_gap(length, coverage)
        if self.candidate.pending_gap_bases > self.merge_distance:
            self.finalize(output)

    def finalize(self, output: List[dict]) -> None:
        candidate = self.candidate
        self.candidate = None
        if candidate is None:
            return

        size = candidate.size
        if size < self.min_size or size > self.max_size:
            return

        coverage = candidate.weighted_coverage / size if size > 0 else 0.0
        ratio = coverage / self.baseline if self.baseline > 0 else math.inf
        low_fraction = candidate.low_bases / size if size > 0 else 0.0

        # Novel candidates are defined by a low-coverage scan. This additional
        # check protects against a merged interval whose normal gap dominates.
        if ratio >= self.low_factor:
            return

        output.append(
            {
                "chrom": candidate.chrom,
                "start": candidate.start,
                "end": candidate.end,
                "size": size,
                "coverage": coverage,
                "chrom_coverage": self.baseline,
                "proportion": ratio,
                "low_fraction": low_fraction,
                "max_low_run": candidate.max_low_run,
            }
        )


@dataclass
class KnownSweep:
    intervals: Sequence[KnownRD]
    next_index: int = 0
    active_indices: List[int] = None  # type: ignore[assignment]
    last_start: Optional[int] = None

    def __post_init__(self) -> None:
        if self.active_indices is None:
            self.active_indices = []


# Initialized once in every worker process.
_WORKER_KNOWN_BY_CHROM: Dict[str, Tuple[KnownRD, ...]] = {}
_WORKER_CONFIG: Dict[str, object] = {}


def init_worker(known_by_chrom: Dict[str, Tuple[KnownRD, ...]], config: Dict[str, object]) -> None:
    global _WORKER_KNOWN_BY_CHROM, _WORKER_CONFIG
    _WORKER_KNOWN_BY_CHROM = known_by_chrom
    _WORKER_CONFIG = config


def open_coverage(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="strict")
    return open(path, "rt", encoding="utf-8", errors="strict")


def iter_coverage(path: str) -> Iterator[Tuple[str, int, int, float]]:
    with open_coverage(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: expected at least 4 tab-separated columns")
            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
                coverage = float(fields[3])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: invalid start/end/coverage: {raw.rstrip()}") from exc
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid interval {chrom}:{start}-{end}")
            if not math.isfinite(coverage) or coverage < 0:
                raise ValueError(f"{path}:{line_number}: invalid coverage {coverage}")
            yield chrom, start, end, coverage


def sample_name_from_path(path: str) -> str:
    name = Path(path).name
    for suffix in (".per-base.bed.gz", ".per-base.bed"):
        if name.endswith(suffix):
            sample = name[: -len(suffix)]
            if sample:
                return sample
    raise ValueError(
        f"Cannot derive sample name from {name!r}; expected *.per-base.bed or *.per-base.bed.gz"
    )


def first_pass(path: str) -> Tuple[Dict[str, ChromStats], int]:
    stats: Dict[str, ChromStats] = {}
    total_records = 0
    previous_end: Dict[str, int] = {}

    for chrom, start, end, coverage in iter_coverage(path):
        prev = previous_end.get(chrom)
        if prev is not None and start < prev:
            raise ValueError(
                f"Coverage BED is not sorted or contains overlapping intervals on {chrom}: "
                f"previous end={prev}, next start={start}"
            )
        previous_end[chrom] = end
        stats.setdefault(chrom, ChromStats()).add(start, end, coverage)
        total_records += 1

    if total_records == 0:
        raise ValueError("Coverage file is empty")
    return stats, total_records


def update_known_intervals(
    chrom: str,
    start: int,
    end: int,
    coverage: float,
    low_threshold: float,
    sweep: KnownSweep,
    accumulators: List[KnownAccumulator],
) -> None:
    if sweep.last_start is not None and start < sweep.last_start:
        raise ValueError(f"Second pass is not sorted on {chrom}: {start} < {sweep.last_start}")
    sweep.last_start = start

    if sweep.active_indices:
        sweep.active_indices = [
            idx for idx in sweep.active_indices if sweep.intervals[idx].end > start
        ]

    while (
        sweep.next_index < len(sweep.intervals)
        and sweep.intervals[sweep.next_index].start < end
    ):
        rd = sweep.intervals[sweep.next_index]
        if rd.end > start:
            sweep.active_indices.append(sweep.next_index)
        sweep.next_index += 1

    for idx in sweep.active_indices:
        rd = sweep.intervals[idx]
        overlap_start = max(start, rd.start)
        overlap_end = min(end, rd.end)
        if overlap_end > overlap_start:
            accumulators[idx].add(overlap_start, overlap_end, coverage, low_threshold)


def classify_ratio(ratio: float, config: Dict[str, object]) -> Optional[str]:
    if ratio <= float(config["complete_threshold"]):
        return "Complete deletion"
    if ratio <= float(config["strong_threshold"]):
        return "Strong deletion"
    if ratio <= float(config["moderate_threshold"]):
        return "Moderate deletion"
    if ratio <= float(config["partial_threshold"]):
        return "Partial deletion"
    return None


def overlap_length(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def process_sample(path: str) -> dict:
    started = time.monotonic()
    sample = sample_name_from_path(path)
    config = _WORKER_CONFIG
    known_by_chrom = _WORKER_KNOWN_BY_CHROM

    chrom_stats, total_records = first_pass(path)
    baselines = {chrom: stat.mean for chrom, stat in chrom_stats.items()}

    primary_chrom = str(config["primary_chrom"])
    if primary_chrom not in baselines:
        raise ValueError(f"Primary chromosome {primary_chrom!r} is absent from coverage file")
    if baselines[primary_chrom] < float(config["min_primary_coverage"]):
        raise ValueError(
            f"Primary chromosome mean coverage {baselines[primary_chrom]:.4f} is below "
            f"--min-primary-coverage={config['min_primary_coverage']}"
        )

    # Allocate only for contigs actually present in the mosdepth file. This
    # prevents absent RvD1/TbD1 contigs from being reported as deletions.
    sweeps: Dict[str, KnownSweep] = {}
    known_accums: Dict[str, List[KnownAccumulator]] = {}
    novel_detectors: Dict[str, NovelDetector] = {}

    for chrom, baseline in baselines.items():
        intervals = known_by_chrom.get(chrom, ())
        if intervals and baseline > 0:
            sweeps[chrom] = KnownSweep(intervals=intervals)
            known_accums[chrom] = [KnownAccumulator() for _ in intervals]
        if baseline > 0:
            novel_detectors[chrom] = NovelDetector(
                chrom=chrom,
                baseline=baseline,
                low_factor=float(config["novel_threshold"]),
                merge_distance=int(config["merge_distance"]),
                min_size=int(config["min_size"]),
                max_size=int(config["max_size"]),
            )

    novel_candidates: List[dict] = []
    last_end_by_chrom: Dict[str, int] = {}

    def handle_segment(chrom: str, start: int, end: int, coverage: float) -> None:
        baseline = baselines.get(chrom, 0.0)
        if baseline <= 0:
            return
        low_threshold = baseline * float(config["novel_threshold"])

        sweep = sweeps.get(chrom)
        if sweep is not None:
            update_known_intervals(
                chrom,
                start,
                end,
                coverage,
                low_threshold,
                sweep,
                known_accums[chrom],
            )

        detector = novel_detectors.get(chrom)
        if detector is not None:
            detector.feed(start, end, coverage, novel_candidates)

    for chrom, start, end, coverage in iter_coverage(path):
        previous_end = last_end_by_chrom.get(chrom)
        if previous_end is not None:
            if start < previous_end:
                raise ValueError(
                    f"Coverage BED is not sorted or overlaps on {chrom}: {start} < {previous_end}"
                )
            if start > previous_end:
                # Mosdepth normally emits zero-coverage blocks. If a gap is
                # nevertheless absent from the file, treat it as zero coverage.
                handle_segment(chrom, previous_end, start, 0.0)
        handle_segment(chrom, start, end, coverage)
        last_end_by_chrom[chrom] = end

    for detector in novel_detectors.values():
        detector.finalize(novel_candidates)

    known_metrics: Dict[Tuple[str, str, int, int], dict] = {}
    known_rows: Dict[Tuple[str, str, int, int], dict] = {}

    for chrom, intervals in known_by_chrom.items():
        if chrom not in chrom_stats or chrom not in known_accums:
            continue
        baseline = baselines[chrom]
        if baseline <= 0:
            continue

        for idx, rd in enumerate(intervals):
            acc = known_accums[chrom][idx]
            size = rd.size
            if size <= 0:
                continue

            missing_bases = max(0, size - acc.covered_bases)
            coverage = acc.weighted_coverage / size
            ratio = coverage / baseline
            low_bases = min(size, acc.low_bases + missing_bases)
            low_fraction = low_bases / size
            max_low_run = max(acc.max_low_run, missing_bases)
            status = classify_ratio(ratio, config)

            key = (chrom, rd.name, rd.start, rd.end)
            metric = {
                "rd": rd,
                "coverage": coverage,
                "chrom_coverage": baseline,
                "proportion": ratio,
                "low_fraction": low_fraction,
                "max_low_run": max_low_run,
                "status": status,
                "detected_start": None,
                "detected_end": None,
            }
            known_metrics[key] = metric

            if status is not None:
                known_rows[key] = {
                    "Sample": sample,
                    "Chromosome": chrom,
                    "Start": rd.start,
                    "End": rd.end,
                    "Size": size,
                    "Type": "DEL",
                    "Source": "known",
                    "KnownRD": rd.name,
                    "ReferenceStart": rd.start,
                    "ReferenceEnd": rd.end,
                    "DetectedStart": None,
                    "DetectedEnd": None,
                    "Coverage": coverage,
                    "ChromCoverage": baseline,
                    "Proportion": ratio,
                    "LowCoverageFraction": low_fraction,
                    "MaxLowCoverageRun": max_low_run,
                    "Status": status,
                    "DetectionReason": "known_interval_ratio",
                    "InputFile": os.path.abspath(path),
                }

    novel_rows: List[dict] = []
    known_overlap_threshold = float(config["known_overlap_threshold"])

    for candidate in novel_candidates:
        best_key: Optional[Tuple[str, str, int, int]] = None
        best_score = 0.0
        chrom = candidate["chrom"]

        for rd in known_by_chrom.get(chrom, ()):
            overlap = overlap_length(candidate["start"], candidate["end"], rd.start, rd.end)
            if overlap <= 0:
                continue
            # Reciprocal-overlap score: both intervals must agree reasonably
            # well; a small RD inside a very broad candidate does not score 1.
            score = overlap / max(candidate["size"], rd.size)
            if score > best_score:
                best_score = score
                best_key = (chrom, rd.name, rd.start, rd.end)

        if best_key is not None and best_score >= known_overlap_threshold:
            metric = known_metrics.get(best_key)
            if metric is not None:
                metric["detected_start"] = candidate["start"]
                metric["detected_end"] = candidate["end"]

                row = known_rows.get(best_key)
                if row is None:
                    # A concentrated low-coverage subinterval can indicate a
                    # partial known deletion even when the mean over the full
                    # expected RD is above --partial-threshold.
                    rd = metric["rd"]
                    row = {
                        "Sample": sample,
                        "Chromosome": chrom,
                        "Start": rd.start,
                        "End": rd.end,
                        "Size": rd.size,
                        "Type": "DEL",
                        "Source": "known",
                        "KnownRD": rd.name,
                        "ReferenceStart": rd.start,
                        "ReferenceEnd": rd.end,
                        "DetectedStart": candidate["start"],
                        "DetectedEnd": candidate["end"],
                        "Coverage": metric["coverage"],
                        "ChromCoverage": metric["chrom_coverage"],
                        "Proportion": metric["proportion"],
                        "LowCoverageFraction": metric["low_fraction"],
                        "MaxLowCoverageRun": metric["max_low_run"],
                        "Status": "Partial deletion",
                        "DetectionReason": "known_overlap_low_coverage_scan",
                        "InputFile": os.path.abspath(path),
                    }
                    known_rows[best_key] = row
                else:
                    row["DetectedStart"] = candidate["start"]
                    row["DetectedEnd"] = candidate["end"]
            continue

        status = classify_ratio(candidate["proportion"], config) or "Moderate deletion"
        novel_rows.append(
            {
                "Sample": sample,
                "Chromosome": chrom,
                "Start": candidate["start"],
                "End": candidate["end"],
                "Size": candidate["size"],
                "Type": "DEL",
                "Source": "novel",
                "KnownRD": "",
                "ReferenceStart": None,
                "ReferenceEnd": None,
                "DetectedStart": candidate["start"],
                "DetectedEnd": candidate["end"],
                "Coverage": candidate["coverage"],
                "ChromCoverage": candidate["chrom_coverage"],
                "Proportion": candidate["proportion"],
                "LowCoverageFraction": candidate["low_fraction"],
                "MaxLowCoverageRun": candidate["max_low_run"],
                "Status": status,
                "DetectionReason": "low_coverage_scan",
                "InputFile": os.path.abspath(path),
            }
        )

    rows = list(known_rows.values()) + novel_rows
    rows.sort(key=lambda row: (row["Chromosome"], int(row["Start"]), int(row["End"]), row["KnownRD"]))

    return {
        "sample": sample,
        "input_file": os.path.abspath(path),
        "rows": rows,
        "records": total_records,
        "primary_coverage": baselines[primary_chrom],
        "elapsed": time.monotonic() - started,
        "novel_candidates": len(novel_candidates),
    }


def load_known_rds(path: str) -> Dict[str, Tuple[KnownRD, ...]]:
    grouped: Dict[str, List[KnownRD]] = {}
    seen: set[Tuple[str, int, int, str]] = set()

    with open(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: expected chrom, start, end, name")
            chrom, start_text, end_text, name = fields[:4]
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError as exc:
                # Permit one conventional header row.
                if line_number == 1 and start_text.lower() == "start":
                    continue
                raise ValueError(f"{path}:{line_number}: invalid coordinates") from exc
            if end <= start:
                raise ValueError(f"{path}:{line_number}: end must be greater than start")
            if name == "H37Rv":
                continue
            key = (chrom, start, end, name)
            if key in seen:
                continue
            seen.add(key)
            grouped.setdefault(chrom, []).append(KnownRD(chrom, start, end, name))

    if not grouped:
        raise ValueError(f"No known RD intervals loaded from {path}")

    return {
        chrom: tuple(sorted(intervals, key=lambda rd: (rd.start, rd.end, rd.name)))
        for chrom, intervals in grouped.items()
    }


def is_coverage_filename(name: str) -> bool:
    return name.endswith(".per-base.bed.gz") or name.endswith(".per-base.bed")


def iter_input_files(input_path: str, recursive: bool) -> Iterator[str]:
    path = os.path.abspath(input_path)
    if os.path.isfile(path):
        if not is_coverage_filename(os.path.basename(path)):
            raise ValueError("Input file must end with .per-base.bed or .per-base.bed.gz")
        yield path
        return
    if not os.path.isdir(path):
        raise FileNotFoundError(path)

    if not recursive:
        with os.scandir(path) as entries:
            for entry in entries:
                if entry.is_file(follow_symlinks=True) and is_coverage_filename(entry.name):
                    yield entry.path
        return

    stack = [path]
    while stack:
        current = stack.pop()
        with os.scandir(current) as entries:
            for entry in entries:
                if entry.is_dir(follow_symlinks=False):
                    stack.append(entry.path)
                elif entry.is_file(follow_symlinks=True) and is_coverage_filename(entry.name):
                    yield entry.path


def count_input_files(input_path: str, recursive: bool) -> int:
    return sum(1 for _ in iter_input_files(input_path, recursive))


def count_input_status(
    input_path: str, recursive: bool, done_samples: set[str]
) -> Tuple[int, int, int]:
    """Return (total, already_done, pending) in one directory scan."""
    total = 0
    already_done = 0
    pending = 0
    for path in iter_input_files(input_path, recursive):
        total += 1
        if sample_name_from_path(path) in done_samples:
            already_done += 1
        else:
            pending += 1
    return total, already_done, pending


def format_duration(seconds: float) -> str:
    if not math.isfinite(seconds) or seconds < 0:
        return "--"
    total = int(round(seconds))
    days, remainder = divmod(total, 86400)
    hours, remainder = divmod(remainder, 3600)
    minutes, secs = divmod(remainder, 60)
    if days:
        return f"{days}d {hours:02d}:{minutes:02d}:{secs:02d}"
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


class ProgressReporter:
    """TTY progress bar with a log-friendly fallback and rolling speed/ETA."""

    def __init__(
        self,
        total: int,
        total_files: int,
        already_done: int,
        interval: float,
        enabled: bool,
        rolling_window: float = 60.0,
        width: int = 28,
    ) -> None:
        self.total = max(0, total)
        self.total_files = max(0, total_files)
        self.already_done = max(0, already_done)
        self.interval = interval
        self.enabled = enabled
        self.rolling_window = rolling_window
        self.width = width
        self.started = time.monotonic()
        self.last_report = 0.0
        self.last_processed = -1
        self.history: deque[Tuple[float, int]] = deque()
        self.is_tty = bool(enabled and sys.stderr.isatty())
        self.line_visible = False

    def clear_line(self) -> None:
        if self.is_tty and self.line_visible:
            print("\r\033[K", end="", file=sys.stderr, flush=True)
            self.line_visible = False

    def _rates(self, now: float, processed: int) -> Tuple[float, float]:
        elapsed = max(0.0, now - self.started)
        average = processed / elapsed if elapsed > 0 else 0.0

        self.history.append((now, processed))
        cutoff = now - self.rolling_window
        while len(self.history) > 1 and self.history[1][0] <= cutoff:
            self.history.popleft()

        recent = 0.0
        if len(self.history) >= 2:
            first_time, first_processed = self.history[0]
            span = now - first_time
            if span > 0:
                recent = (processed - first_processed) / span
        return max(0.0, recent), max(0.0, average)

    def update(
        self,
        processed: int,
        done: int,
        errors: int,
        deletion_rows: int,
        queued: int,
        force: bool = False,
    ) -> None:
        if not self.enabled:
            return

        now = time.monotonic()
        terminal_state = self.total == 0 or processed >= self.total
        if (
            not force
            and processed != 1
            and not terminal_state
            and now - self.last_report < self.interval
        ):
            return

        recent_rate, average_rate = self._rates(now, processed)
        elapsed = now - self.started
        remaining = max(0, self.total - processed)

        # Use rolling speed after enough observations; fall back to run average
        # during startup. If the rolling window has stalled, ETA is unknown.
        history_span = self.history[-1][0] - self.history[0][0] if len(self.history) > 1 else 0.0
        if history_span >= min(10.0, self.rolling_window / 2):
            eta_rate = recent_rate
        else:
            eta_rate = average_rate
        eta = remaining / eta_rate if eta_rate > 0 else math.inf

        fraction = processed / self.total if self.total > 0 else 1.0
        fraction = min(1.0, max(0.0, fraction))
        filled = int(round(self.width * fraction))
        bar = "#" * filled + "-" * (self.width - filled)
        percent = 100.0 * fraction
        overall_done = min(self.total_files, self.already_done + done)

        line = (
            f"Progress [{bar}] {percent:6.2f}%  {processed}/{self.total} "
            f"| ok={done} err={errors} rows={deletion_rows} "
            f"| speed={recent_rate:.2f}/s avg={average_rate:.2f}/s "
            f"| ETA={format_duration(eta)} elapsed={format_duration(elapsed)} "
            f"| queued={queued} overall_done={overall_done}/{self.total_files}"
        )

        if self.is_tty:
            print(f"\r{line}\033[K", end="", file=sys.stderr, flush=True)
            self.line_visible = True
        else:
            LOG.info("%s", line)

        self.last_report = now
        self.last_processed = processed

    def finish(
        self, processed: int, done: int, errors: int, deletion_rows: int, queued: int = 0
    ) -> None:
        # Avoid duplicate final lines in redirected logs. In an interactive
        # terminal, still terminate the in-place progress line with a newline.
        if self.last_processed != processed:
            self.update(processed, done, errors, deletion_rows, queued, force=True)
        if self.is_tty and self.line_visible:
            print(file=sys.stderr, flush=True)
            self.line_visible = False


OUTPUT_COLUMNS = [
    "Sample",
    "Chromosome",
    "Start",
    "End",
    "Size",
    "Type",
    "Source",
    "KnownRD",
    "ReferenceStart",
    "ReferenceEnd",
    "DetectedStart",
    "DetectedEnd",
    "Coverage",
    "ChromCoverage",
    "Proportion",
    "LowCoverageFraction",
    "MaxLowCoverageRun",
    "Status",
    "DetectionReason",
]


def connect_state_db(path: str) -> sqlite3.Connection:
    connection = sqlite3.connect(path, timeout=60.0)
    connection.execute("PRAGMA journal_mode=WAL")
    connection.execute("PRAGMA synchronous=NORMAL")
    connection.execute("PRAGMA temp_store=MEMORY")
    connection.execute(
        """
        CREATE TABLE IF NOT EXISTS metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        )
        """
    )
    connection.execute(
        """
        CREATE TABLE IF NOT EXISTS samples (
            sample TEXT PRIMARY KEY,
            input_file TEXT NOT NULL,
            status TEXT NOT NULL,
            error TEXT,
            records INTEGER,
            primary_coverage REAL,
            elapsed REAL,
            processed_at TEXT DEFAULT CURRENT_TIMESTAMP
        )
        """
    )
    connection.execute(
        """
        CREATE TABLE IF NOT EXISTS deletions (
            sample TEXT NOT NULL,
            chromosome TEXT NOT NULL,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            size INTEGER NOT NULL,
            type TEXT NOT NULL,
            source TEXT NOT NULL,
            known_rd TEXT NOT NULL DEFAULT '',
            reference_start INTEGER,
            reference_end INTEGER,
            detected_start INTEGER,
            detected_end INTEGER,
            coverage REAL NOT NULL,
            chrom_coverage REAL NOT NULL,
            proportion REAL NOT NULL,
            low_coverage_fraction REAL NOT NULL,
            max_low_coverage_run INTEGER NOT NULL,
            status TEXT NOT NULL,
            detection_reason TEXT NOT NULL,
            input_file TEXT NOT NULL,
            PRIMARY KEY (sample, chromosome, start, end, source, known_rd)
        )
        """
    )
    connection.execute("CREATE INDEX IF NOT EXISTS deletions_sample_idx ON deletions(sample)")
    connection.commit()
    return connection



def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ensure_run_signature(connection: sqlite3.Connection, signature: str) -> None:
    row = connection.execute(
        "SELECT value FROM metadata WHERE key='run_signature'"
    ).fetchone()
    if row is None:
        with connection:
            connection.execute(
                "INSERT INTO metadata(key, value) VALUES('run_signature', ?)",
                (signature,),
            )
        return
    if row[0] != signature:
        raise RuntimeError(
            "Existing resume database was created with different RD definitions or "
            "thresholds. Use --fresh or choose a different --output path."
        )


def completed_samples(connection: sqlite3.Connection) -> set[str]:
    return {row[0] for row in connection.execute("SELECT sample FROM samples WHERE status='done'")}


def save_success(connection: sqlite3.Connection, result: dict) -> None:
    rows = result["rows"]
    with connection:
        connection.execute("DELETE FROM deletions WHERE sample=?", (result["sample"],))
        if rows:
            connection.executemany(
                """
                INSERT OR REPLACE INTO deletions (
                    sample, chromosome, start, end, size, type, source, known_rd,
                    reference_start, reference_end, detected_start, detected_end,
                    coverage, chrom_coverage, proportion, low_coverage_fraction,
                    max_low_coverage_run, status, detection_reason, input_file
                ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
                """,
                [
                    (
                        row["Sample"], row["Chromosome"], row["Start"], row["End"], row["Size"],
                        row["Type"], row["Source"], row["KnownRD"], row["ReferenceStart"],
                        row["ReferenceEnd"], row["DetectedStart"], row["DetectedEnd"],
                        row["Coverage"], row["ChromCoverage"], row["Proportion"],
                        row["LowCoverageFraction"], row["MaxLowCoverageRun"], row["Status"],
                        row["DetectionReason"], row["InputFile"],
                    )
                    for row in rows
                ],
            )
        connection.execute(
            """
            INSERT OR REPLACE INTO samples
                (sample, input_file, status, error, records, primary_coverage, elapsed, processed_at)
            VALUES (?, ?, 'done', NULL, ?, ?, ?, CURRENT_TIMESTAMP)
            """,
            (
                result["sample"], result["input_file"], result["records"],
                result["primary_coverage"], result["elapsed"],
            ),
        )


def save_error(connection: sqlite3.Connection, sample: str, input_file: str, error: str) -> None:
    with connection:
        connection.execute("DELETE FROM deletions WHERE sample=?", (sample,))
        connection.execute(
            """
            INSERT OR REPLACE INTO samples
                (sample, input_file, status, error, records, primary_coverage, elapsed, processed_at)
            VALUES (?, ?, 'error', ?, NULL, NULL, NULL, CURRENT_TIMESTAMP)
            """,
            (sample, os.path.abspath(input_file), error[:10000]),
        )


def format_tsv_value(column: str, value: object) -> object:
    if value is None:
        return ""
    if column in {"Coverage", "ChromCoverage"}:
        return f"{float(value):.6f}"
    if column in {"Proportion", "LowCoverageFraction"}:
        return f"{float(value):.8f}"
    return value


def export_tsv(connection: sqlite3.Connection, output_path: str) -> None:
    temporary = output_path + ".tmp"
    query = """
        SELECT sample, chromosome, start, end, size, type, source, known_rd,
               reference_start, reference_end, detected_start, detected_end,
               coverage, chrom_coverage, proportion, low_coverage_fraction,
               max_low_coverage_run, status, detection_reason
        FROM deletions
        ORDER BY sample, chromosome, start, end, source, known_rd
    """
    with open(temporary, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(OUTPUT_COLUMNS)
        for values in connection.execute(query):
            writer.writerow(
                [format_tsv_value(column, value) for column, value in zip(OUTPUT_COLUMNS, values)]
            )
    os.replace(temporary, output_path)


def export_errors(connection: sqlite3.Connection, path: str) -> None:
    temporary = path + ".tmp"
    with open(temporary, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Sample", "InputFile", "Error"])
        for row in connection.execute(
            "SELECT sample, input_file, error FROM samples WHERE status='error' ORDER BY sample"
        ):
            writer.writerow(row)
    os.replace(temporary, path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Batch known+novel RD detection from mosdepth per-base BED files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i", "--input", required=True, help="One *.per-base.bed[.gz] file or a directory")
    parser.add_argument("-k", "--known-rds", required=True, help="Four-column BED-like file: chrom, start, end, RD name")
    parser.add_argument("-o", "--output", required=True, help="Final combined TSV")
    parser.add_argument("-j", "--jobs", type=int, default=min(8, os.cpu_count() or 1), help="Parallel worker processes")
    parser.add_argument("--max-pending", type=int, default=0, help="Maximum queued samples; 0 means 2 × jobs")
    parser.add_argument("--recursive", action="store_true", help="Search input directory recursively")
    parser.add_argument("--limit", type=int, default=0, help="Process at most this many pending samples; 0 means all")
    parser.add_argument("--fresh", action="store_true", help="Delete previous state and recompute from scratch")
    parser.add_argument(
        "--progress-interval",
        type=float,
        default=10.0,
        help="Seconds between progress updates; terminal output uses one updating line",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress, speed and ETA reporting",
    )

    parser.add_argument("--primary-chrom", default="NC_000962.3", help="Primary chromosome required in every sample")
    parser.add_argument("--min-primary-coverage", type=float, default=1.0, help="Reject samples below this chromosome-wide mean coverage")
    parser.add_argument("--novel-threshold", type=float, default=0.10, help="Novel scan: coverage below this fraction of chromosome mean is low")
    parser.add_argument("--merge-distance", type=int, default=100, help="Maximum normal-coverage gap merged between low-coverage blocks")
    parser.add_argument("--min-size", type=int, default=500, help="Minimum novel deletion size")
    parser.add_argument("--max-size", type=int, default=100000, help="Maximum novel deletion size")
    parser.add_argument("--known-overlap-threshold", type=float, default=0.50, help="Reciprocal overlap needed to annotate a candidate as a known RD")

    parser.add_argument("--complete-threshold", type=float, default=0.01, help="Coverage ratio at or below this is Complete deletion")
    parser.add_argument("--strong-threshold", type=float, default=0.05, help="Coverage ratio at or below this is Strong deletion")
    parser.add_argument("--moderate-threshold", type=float, default=0.10, help="Coverage ratio at or below this is Moderate deletion")
    parser.add_argument("--partial-threshold", type=float, default=0.50, help="Coverage ratio at or below this is Partial deletion")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if args.jobs < 1:
        raise ValueError("--jobs must be at least 1")
    if args.max_pending < 0:
        raise ValueError("--max-pending cannot be negative")
    if args.limit < 0:
        raise ValueError("--limit cannot be negative")
    if args.progress_interval <= 0:
        raise ValueError("--progress-interval must be greater than 0")
    if args.min_size < 1 or args.max_size < args.min_size:
        raise ValueError("Require 1 <= --min-size <= --max-size")
    if args.merge_distance < 0:
        raise ValueError("--merge-distance cannot be negative")
    if args.min_primary_coverage < 0:
        raise ValueError("--min-primary-coverage cannot be negative")
    for name in (
        "novel_threshold", "known_overlap_threshold", "complete_threshold",
        "strong_threshold", "moderate_threshold", "partial_threshold",
    ):
        value = float(getattr(args, name))
        if not 0 <= value <= 1:
            raise ValueError(f"--{name.replace('_', '-')} must be between 0 and 1")
    if not (
        args.complete_threshold
        <= args.strong_threshold
        <= args.moderate_threshold
        <= args.partial_threshold
    ):
        raise ValueError("Deletion thresholds must be nondecreasing: complete <= strong <= moderate <= partial")


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        validate_args(args)
    except ValueError as exc:
        parser.error(str(exc))

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(message)s",
    )

    output_path = os.path.abspath(args.output)
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    state_path = output_path + ".state.sqlite"
    errors_path = output_path + ".errors.tsv"

    if args.fresh:
        for path in (state_path, state_path + "-wal", state_path + "-shm", output_path, errors_path):
            try:
                os.remove(path)
            except FileNotFoundError:
                pass

    try:
        known_by_chrom = load_known_rds(args.known_rds)
    except Exception as exc:
        LOG.error("Initialization failed: %s", exc)
        return 2

    connection = connect_state_db(state_path)

    config: Dict[str, object] = {
        "primary_chrom": args.primary_chrom,
        "min_primary_coverage": args.min_primary_coverage,
        "novel_threshold": args.novel_threshold,
        "merge_distance": args.merge_distance,
        "min_size": args.min_size,
        "max_size": args.max_size,
        "known_overlap_threshold": args.known_overlap_threshold,
        "complete_threshold": args.complete_threshold,
        "strong_threshold": args.strong_threshold,
        "moderate_threshold": args.moderate_threshold,
        "partial_threshold": args.partial_threshold,
    }

    signature_payload = {
        "script_version": 1,
        "known_rds_sha256": sha256_file(args.known_rds),
        "config": config,
    }
    run_signature = hashlib.sha256(
        json.dumps(signature_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    try:
        ensure_run_signature(connection, run_signature)
    except RuntimeError as exc:
        LOG.error("%s", exc)
        connection.close()
        return 2

    done_samples = completed_samples(connection)
    try:
        total_files, already_done_in_input, pending_count = count_input_status(
            args.input, args.recursive, done_samples
        )
    except Exception as exc:
        LOG.error("Input scan failed: %s", exc)
        connection.close()
        return 2

    if total_files == 0:
        LOG.error("No *.per-base.bed or *.per-base.bed.gz files found")
        connection.close()
        return 2


    run_target = min(pending_count, args.limit) if args.limit else pending_count
    max_pending = args.max_pending or args.jobs * 2
    pending_generator = (
        path
        for path in iter_input_files(args.input, args.recursive)
        if sample_name_from_path(path) not in done_samples
    )

    LOG.info(
        "Found %d coverage files; %d already complete in this input; "
        "%d pending; this run target=%d; workers=%d; max_pending=%d",
        total_files,
        already_done_in_input,
        pending_count,
        run_target,
        args.jobs,
        max_pending,
    )

    progress = ProgressReporter(
        total=run_target,
        total_files=total_files,
        already_done=already_done_in_input,
        interval=args.progress_interval,
        enabled=not args.no_progress,
    )

    submitted = 0
    completed_now = 0
    errors_now = 0
    detection_rows_now = 0
    started = time.monotonic()
    futures = {}

    def next_path() -> Optional[str]:
        nonlocal submitted
        if args.limit and submitted >= args.limit:
            return None
        try:
            path = next(pending_generator)
        except StopIteration:
            return None
        submitted += 1
        return path

    try:
        with ProcessPoolExecutor(
            max_workers=args.jobs,
            initializer=init_worker,
            initargs=(known_by_chrom, config),
        ) as executor:
            for _ in range(max_pending):
                path = next_path()
                if path is None:
                    break
                futures[executor.submit(process_sample, path)] = path

            progress.update(
                processed=0,
                done=0,
                errors=0,
                deletion_rows=0,
                queued=len(futures),
                force=True,
            )

            while futures:
                finished, _ = wait(futures, timeout=1.0, return_when=FIRST_COMPLETED)
                if not finished:
                    progress.update(
                        processed=completed_now + errors_now,
                        done=completed_now,
                        errors=errors_now,
                        deletion_rows=detection_rows_now,
                        queued=len(futures),
                    )
                    continue

                for future in finished:
                    path = futures.pop(future)
                    try:
                        result = future.result()
                        save_success(connection, result)
                        completed_now += 1
                        detection_rows_now += len(result["rows"])
                    except Exception as exc:
                        try:
                            sample = sample_name_from_path(path)
                        except Exception:
                            sample = Path(path).name
                        error_text = f"{type(exc).__name__}: {exc}"
                        save_error(connection, sample, path, error_text)
                        errors_now += 1
                        progress.clear_line()
                        LOG.error("%s: %s", sample, error_text)

                    replacement = next_path()
                    if replacement is not None:
                        futures[executor.submit(process_sample, replacement)] = replacement

                progress.update(
                    processed=completed_now + errors_now,
                    done=completed_now,
                    errors=errors_now,
                    deletion_rows=detection_rows_now,
                    queued=len(futures),
                )

            progress.finish(
                processed=completed_now + errors_now,
                done=completed_now,
                errors=errors_now,
                deletion_rows=detection_rows_now,
            )
    except KeyboardInterrupt:
        progress.clear_line()
        LOG.warning("Interrupted. Completed samples are safely stored in %s", state_path)
        export_tsv(connection, output_path)
        export_errors(connection, errors_path)
        connection.close()
        return 130
    except Exception as exc:
        progress.clear_line()
        LOG.exception("Batch execution failed: %s", exc)
        export_tsv(connection, output_path)
        export_errors(connection, errors_path)
        connection.close()
        return 1

    export_tsv(connection, output_path)
    export_errors(connection, errors_path)

    totals = connection.execute(
        "SELECT status, COUNT(*) FROM samples GROUP BY status"
    ).fetchall()
    total_deletions = connection.execute("SELECT COUNT(*) FROM deletions").fetchone()[0]
    connection.close()

    elapsed = time.monotonic() - started
    LOG.info("Finished in %.1f s", elapsed)
    LOG.info("Sample states: %s", dict(totals))
    LOG.info("Final deletion rows: %d", total_deletions)
    LOG.info("Combined TSV: %s", output_path)
    LOG.info("Errors TSV: %s", errors_path)
    LOG.info("Resume state: %s", state_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
