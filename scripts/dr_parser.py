#!/usr/bin/env python3
"""
Сбор признаков лекарственной устойчивости из TB‑Profiler JSON‑отчётов.

Пример использования
--------------------
# вывод в файл
python dr_parser.py -o result.tsv sample1.json sample2.json

# вывод на экран
python dr_parser.py sample.json
"""

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import TextIO, Tuple

COLUMNS = [
    "Sample",
    "Rifampicin",
    "Isoniazid",
    "Ethambutol",
    "Pyrazinamide",
    "Moxifloxacin",
    "Levofloxacin",
    "Bedaquiline",
    "Delamanid",
    "Pretomanid",
    "Linezolid",
    "Streptomycin",
    "Amikacin",
    "Kanamycin",
    "Capreomycin",
    "Clofazimine",
    "Ethionamide",
    "Para-aminosalicylic_acid",
    "Cycloserine",
]

# соответствие «имя из JSON → имя колонки»
JSON2COL = {c.lower(): c for c in COLUMNS[1:]}


def parse_file(path: Path) -> Tuple[str, dict[str, str]]:
    """
    Возвращает (sample_id, row_dict),
    где row_dict: колонка → '' или 'Assoc w R'.
    """
    with path.open() as f:
        data = json.load(f)

    sample_id: str = data.get("id", "UNKNOWN")

    row = {col: "" for col in COLUMNS[1:]}  # без Sample
    for var in data.get("dr_variants", []):
        for drug_info in var.get("drugs", []):
            drug = drug_info.get("drug", "").lower()
            if (
                drug in JSON2COL
                and drug_info.get("confidence", "").strip() == "Assoc w R"
            ):
                row[JSON2COL[drug]] = "Assoc w R"

    return sample_id, row


def write_rows(
    rows: list[list[str]], destination: TextIO, write_header: bool = True
) -> None:
    writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
    if write_header:
        writer.writerow(COLUMNS)
    writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parse TB‑Profiler JSON and output TSV with resistance flags."
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="FILE",
        help="Файл для вывода (TSV). Если не указан, вывод в stdout",
    )
    parser.add_argument("json_files", nargs="+", help="Файлы *.json")
    args = parser.parse_args()

    rows = []
    for fp in args.json_files:
        sample_id, parsed = parse_file(Path(fp))
        rows.append([sample_id] + [parsed[col] for col in COLUMNS[1:]])

    if args.output:
        with open(args.output, "w", newline="") as f:
            write_rows(rows, f)
    else:
        write_rows(rows, sys.stdout)


if __name__ == "__main__":
    main()
