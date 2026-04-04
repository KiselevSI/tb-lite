#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import os
from pathlib import Path


REQUIRED_COLS = {
    "Chromosome",
    "Start",
    "End",
    "Size",
    "Type",
    "Coverage",
    "ChromCoverage",
    "Proportion",
    "DeletionCategory",
    "KnownRD",
}

OUT_COLS = [
    "sample_id",
    "start_pos",
    "end_pos",
    "size_bp",
    "type",
    "coverage",
    "chrom_coverage",
    "proportion",
    "known_rd",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Сбор TSV-файлов с делециями (novel_rd.tsv) в один CSV (без БД). "
                    "На вход можно давать директории: будут обработаны рекурсивно."
    )
    parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        help="Файлы и/или директории. Директории обходятся рекурсивно.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Путь к выходному CSV файлу",
    )
    parser.add_argument(
        "--pattern",
        default="*.novel_rd.tsv",
        help="Шаблон файлов для поиска внутри директорий (по умолчанию: *.novel_rd.tsv)",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Не печатать прогресс и предупреждения (только ошибки Python при падении).",
    )
    parser.add_argument(
        "--summary",
        action="store_true",
        help="Печатать только итоговую сводку в конце (stderr).",
    )
    return parser.parse_args()


def log(msg: str, quiet: bool, stream=None):
    if quiet:
        return
    if stream is None:
        stream = os.sys.stderr
    print(msg, file=stream)


def filename_to_sample_id(path: str) -> str:
    return Path(path).name.split(".")[0]


def to_int(value):
    value = (value or "").strip()
    if value == "" or value.lower() == "na":
        return None
    try:
        return int(value)
    except ValueError:
        return None


def to_float(value):
    value = (value or "").strip()
    if value == "" or value.lower() == "na":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def iter_input_files(inputs, pattern: str, quiet: bool):
    for item in inputs:
        p = Path(item)

        if p.is_file():
            yield p
            continue

        if p.is_dir():
            yield from p.rglob(pattern)
            continue

        log(f"[WARN] Не найдено (не файл и не директория): {item}", quiet)


def iter_deletions_from_file(filepath: Path, quiet: bool):
    sample_id = filename_to_sample_id(str(filepath))
    # УБРАН вывод по каждому файлу

    try:
        with filepath.open("r", encoding="utf-8", errors="replace") as f:
            reader = csv.DictReader(f, delimiter="\t")

            if not reader.fieldnames:
                log(f"[WARN] {filepath}: пустой/без заголовка. Пропуск.", quiet)
                return

            missing = REQUIRED_COLS - set(reader.fieldnames)
            if missing:
                log(
                    f"[WARN] {filepath}: нет обязательных колонок {sorted(missing)}. Пропуск.",
                    quiet,
                )
                return

            for line in reader:
                start_pos = to_int(line.get("Start"))
                end_pos = to_int(line.get("End"))
                size_bp = to_int(line.get("Size"))
                type_ = (line.get("Type") or "").strip()

                coverage = to_float(line.get("Coverage"))
                chrom_coverage = to_float(line.get("ChromCoverage"))
                proportion = to_float(line.get("Proportion"))
                known_rd = (line.get("KnownRD") or "").strip() or ""

                if start_pos is None or end_pos is None or size_bp is None or not type_:
                    continue

                yield {
                    "sample_id": sample_id,
                    "start_pos": start_pos,
                    "end_pos": end_pos,
                    "size_bp": size_bp,
                    "type": type_,
                    "coverage": coverage,
                    "chrom_coverage": chrom_coverage,
                    "proportion": proportion,
                    "known_rd": known_rd,
                }

    except Exception as e:
        log(f"[WARN] Ошибка чтения {filepath}: {e}", quiet)


def main():
    args = parse_args()

    total_files = 0
    total_rows = 0

    with open(args.output, "w", encoding="utf-8", newline="") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=OUT_COLS)
        writer.writeheader()

        for fp in iter_input_files(args.input, args.pattern, args.quiet):
            if not fp.exists() or not fp.is_file():
                continue

            total_files += 1
            for row in iter_deletions_from_file(fp, args.quiet):
                writer.writerow(row)
                total_rows += 1

    if args.summary and not args.quiet:
        print(f"[DONE] Файлов обработано: {total_files}, строк записано: {total_rows}.", file=os.sys.stderr)
        print(f"[DONE] CSV: {args.output}", file=os.sys.stderr)


if __name__ == "__main__":
    main()
