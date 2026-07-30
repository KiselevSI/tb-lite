#!/usr/bin/env python3
"""Полные spoligo-таблицы для TB Platform.

SpoTyping запускается с ``--noQuery``, поэтому сам по себе даёт только бинарный
и восьмеричный коды. Этот скрипт добирает всё остальное, что показывает сайт:

* количество ридов на каждый из 43 спейсеров -- из ``<sample>.log``;
* пороги ``min``/``rmin`` -- оттуда же;
* SIT/клада/география -- джойном восьмеричного кода со справочником SpolDB4.

Формат лога SpoTyping::

    ## <input files>
    ## min=5 rmin=6
    ## Spacer	Error-free_number	1-error-tolerant_number	Code
    Spacer1	92	95	1
    ...
    Spacer43	0	0	0

Выходные файлы:

* ``--spacers-output`` -- ``sample_id  spacer1 … spacer43``
  (таблица сайта ``general_spoligo_spacers``);
* ``--full-output``    -- ``Sample  SpolBin  Spol8  MinReads  RminReads
  SIT  Clade  Total  Geography``.

``spotyping.total.tsv`` намеренно остаётся трёхколоночным: его формат ждут
существующие импортёры TB Platform.
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

N_SPACERS = 43

# Индекс колонки в строке лога: 0 = Spacer, 1 = Error-free, 2 = 1-error-tolerant, 3 = Code
COUNT_COLUMN = {"error_free": 1, "tolerant": 2}

MINMAX_RE = re.compile(r"min\s*=\s*(\d+).*?rmin\s*=\s*(\d+)", re.IGNORECASE)
SPACER_ROW_RE = re.compile(r"^Spacer(\d+)\b", re.IGNORECASE)

FULL_COLUMNS = [
    "Sample", "SpolBin", "Spol8", "MinReads", "RminReads",
    "SIT", "Clade", "Total", "Geography",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build the full spoligotyping tables for TB Platform.",
    )
    parser.add_argument(
        "--spoligo-list",
        required=True,
        help="Text file with one per-sample SpoTyping TSV path per line.",
    )
    parser.add_argument(
        "--log-list",
        required=True,
        help="Text file with one SpoTyping .log path per line.",
    )
    parser.add_argument(
        "--spoldb4",
        required=True,
        help="SpolDB4 reference TSV: sit, octal, total, clade, geography.",
    )
    parser.add_argument(
        "--count",
        choices=sorted(COUNT_COLUMN),
        default="tolerant",
        help="Which read count to store per spacer.",
    )
    parser.add_argument(
        "--spacers-output",
        required=True,
        help="Output TSV with per-spacer read counts.",
    )
    parser.add_argument(
        "--full-output",
        required=True,
        help="Output TSV with codes, thresholds and SpolDB4 annotation.",
    )
    return parser.parse_args()


@dataclass
class LogRecord:
    """Разобранный лог SpoTyping одного образца."""

    min_reads: str = ""
    rmin_reads: str = ""
    counts: List[str] = field(default_factory=list)


def read_input_list(path: Path) -> List[Path]:
    items: List[Path] = []
    seen: set[str] = set()

    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            key = str(Path(line))
            if key in seen:
                continue

            seen.add(key)
            items.append(Path(line))

    return items


def sample_id_from_log(path: Path) -> str:
    """ERR13014410.log -> ERR13014410."""
    return path.name[: -len(".log")] if path.name.endswith(".log") else path.name


def parse_log(path: Path, count_idx: int) -> Optional[LogRecord]:
    """Вернуть разбор лога или None, если строк со спейсерами нет."""
    min_reads: Optional[str] = None
    rmin_reads: Optional[str] = None
    counts: Dict[int, str] = {}

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                match = MINMAX_RE.search(line)
                if match:
                    min_reads, rmin_reads = match.group(1), match.group(2)
                continue
            row = SPACER_ROW_RE.match(line)
            if not row:
                continue
            parts = line.split()
            if len(parts) <= count_idx:
                continue
            counts[int(row.group(1))] = parts[count_idx]

    if not counts:
        return None

    return LogRecord(
        min_reads=min_reads or "",
        rmin_reads=rmin_reads or "",
        counts=[counts.get(index, "") for index in range(1, N_SPACERS + 1)],
    )


def read_spoldb4(path: Path) -> Dict[str, Dict[str, str]]:
    """Справочник SpolDB4, ключ -- восьмеричный код."""
    reference: Dict[str, Dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            octal = (row.get("octal") or "").strip()
            if not octal:
                continue
            reference.setdefault(octal, {
                "sit": (row.get("sit") or "").strip(),
                "clade": (row.get("clade") or "").strip(),
                "total": (row.get("total") or "").strip(),
                "geography": (row.get("geography") or "").strip(),
            })
    return reference


def read_codes(paths: List[Path]) -> Dict[str, Dict[str, str]]:
    """Прочитать per-sample таблицы SpoTyping (без шапки): Sample, SpolBin, Spol8."""
    codes: Dict[str, Dict[str, str]] = {}
    for path in paths:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if len(row) < 3 or not row[0].strip():
                    continue
                codes[row[0].strip()] = {
                    "spol_bin": row[1].strip(),
                    "spol8": row[2].strip(),
                }
    return codes


def main() -> None:
    args = parse_args()
    count_idx = COUNT_COLUMN[args.count]

    codes = read_codes(read_input_list(Path(args.spoligo_list)))
    spoldb4 = read_spoldb4(Path(args.spoldb4))

    logs: Dict[str, LogRecord] = {}
    for path in read_input_list(Path(args.log_list)):
        parsed = parse_log(path, count_idx)
        if parsed is not None:
            logs[sample_id_from_log(path)] = parsed

    spacers_path = Path(args.spacers_output)
    with spacers_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id"] + [f"spacer{i}" for i in range(1, N_SPACERS + 1)])
        for sample_id in sorted(logs):
            writer.writerow([sample_id] + logs[sample_id].counts)

    full_path = Path(args.full_output)
    with full_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(FULL_COLUMNS)
        for sample_id in sorted(set(codes) | set(logs)):
            code = codes.get(sample_id, {})
            log = logs.get(sample_id, LogRecord())
            annotation = spoldb4.get(code.get("spol8", ""), {})
            writer.writerow([
                sample_id,
                code.get("spol_bin", ""),
                code.get("spol8", ""),
                log.min_reads,
                log.rmin_reads,
                annotation.get("sit", ""),
                annotation.get("clade", ""),
                annotation.get("total", ""),
                annotation.get("geography", ""),
            ])


if __name__ == "__main__":
    main()
