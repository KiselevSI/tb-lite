#!/usr/bin/env python3
"""Оставить в TSV только строки образцов из списка.

Нужно потому, что метрики покрытия считаются до фильтра качества, и в
general.tsv / tbmix.total.tsv попадают образцы, которые не дошли до вызова
вариантов. В таблицах для TB Platform их быть не должно: в базе они дали бы
карточки без линии, сполиготипа, устойчивости и SNP.

Пример:
    filter_table_by_samples.py -i general.tsv -s passed_samples.txt \\
        --id-column ID -o general.filtered.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def read_samples(path: Path) -> set:
    with open(path, "r", encoding="utf-8") as handle:
        return {line.strip() for line in handle if line.strip()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", required=True, type=Path)
    parser.add_argument("-s", "--samples", required=True, type=Path,
                        help="Файл со списком sample_id, по одному на строку.")
    parser.add_argument("--id-column", required=True,
                        help="Имя колонки с идентификатором образца.")
    parser.add_argument("-o", "--output", required=True, type=Path)
    parser.add_argument("--delimiter", default="\t")
    args = parser.parse_args()

    keep = read_samples(args.samples)

    with open(args.input, "r", encoding="utf-8", newline="") as src:
        reader = csv.reader(src, delimiter=args.delimiter)
        try:
            header = next(reader)
        except StopIteration:
            print(f"ERROR: пустой файл {args.input}", file=sys.stderr)
            return 1

        if args.id_column not in header:
            print(
                f"ERROR: в {args.input} нет колонки {args.id_column!r}; "
                f"есть: {', '.join(header)}",
                file=sys.stderr,
            )
            return 1
        id_index = header.index(args.id_column)

        kept = 0
        dropped = 0
        with open(args.output, "w", encoding="utf-8", newline="") as dst:
            writer = csv.writer(dst, delimiter=args.delimiter, lineterminator="\n")
            writer.writerow(header)
            for row in reader:
                if not row:
                    continue
                if len(row) > id_index and row[id_index].strip() in keep:
                    writer.writerow(row)
                    kept += 1
                else:
                    dropped += 1

    print(
        f"{args.input.name}: оставлено {kept}, отброшено {dropped} "
        f"(список: {len(keep)} образцов)",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
