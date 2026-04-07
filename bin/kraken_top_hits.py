#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path


SAMPLE_FRAC_SUFFIX = ".bracken.add_Unclassified.tsv_frac"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build per-sample top-N organism summary from Kraken/Bracken combined outputs."
    )
    parser.add_argument(
        "-i",
        "--input",
        nargs="*",
        default=[],
        help="Kraken combined *.all_samples.txt files.",
    )
    parser.add_argument(
        "--input-list",
        help="Text file with one Kraken combined path per line.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="kraken.top_hits.tsv",
        help="Output TSV path.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=5,
        help="Number of top organisms per database to include.",
    )
    args = parser.parse_args()

    if not args.input and not args.input_list:
        parser.error("provide --input and/or --input-list")

    return args


def read_input_list(path: str | None) -> list[str]:
    if not path:
        return []

    items: list[str] = []
    with Path(path).open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            items.append(line)

    return items


def merge_inputs(cli_items: list[str], list_path: str | None) -> list[Path]:
    merged: list[Path] = []
    seen: set[str] = set()

    for item in [*(cli_items or []), *read_input_list(list_path)]:
        path = Path(item)
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        merged.append(path)

    return sorted(merged)


def sanitize_label(value: str) -> str:
    label = re.sub(r"[^0-9A-Za-z_]+", "_", value.strip())
    label = re.sub(r"_+", "_", label).strip("_")
    return label or "kraken"


def derive_db_label(path: Path) -> str:
    name = path.name
    if name.endswith(".all_samples.txt"):
        name = name[: -len(".all_samples.txt")]
    else:
        name = path.stem
    return sanitize_label(name)


def sample_from_frac_column(column_name: str) -> str | None:
    if not column_name.endswith(SAMPLE_FRAC_SUFFIX):
        return None
    return column_name[: -len(SAMPLE_FRAC_SUFFIX)] or None


def parse_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("-inf")


def format_top_hit(name: str, frac: str) -> str:
    if not frac:
        return name
    try:
        frac_text = f"{float(frac):.2f}"
    except (TypeError, ValueError):
        frac_text = str(frac)
    return f"{name} ({frac_text})"


def top_hits_for_file(path: Path, top_n: int) -> tuple[list[str], dict[str, list[str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    frac_columns = [name for name in (reader.fieldnames or []) if sample_from_frac_column(name)]
    sample_ids = [sample_from_frac_column(name) for name in frac_columns]

    db_label = derive_db_label(path)
    column_names = [f"{db_label}_top{i}" for i in range(1, top_n + 1)]
    data: dict[str, list[str]] = {}

    for frac_column, sample_id in zip(frac_columns, sample_ids):
        ranked = sorted(
            (
                (parse_float(row.get(frac_column, "")), row.get("name", ""), row.get(frac_column, ""))
                for row in rows
                if row.get("name", "") and row.get(frac_column, "") != ""
            ),
            key=lambda item: item[0],
            reverse=True,
        )
        top_names = [format_top_hit(name, frac) for _score, name, frac in ranked[:top_n]]
        top_names.extend([""] * max(0, top_n - len(top_names)))
        data[sample_id] = top_names

    return column_names, data


def main():
    args = parse_args()

    inputs = merge_inputs(args.input, args.input_list)
    sample_rows: dict[str, dict[str, str]] = {}
    ordered_columns: list[str] = []

    for path in inputs:
        columns, data = top_hits_for_file(path, args.top_n)
        ordered_columns.extend(columns)

        for sample_id, top_names in data.items():
            row = sample_rows.setdefault(sample_id, {"ID": sample_id})
            for column_name, organism in zip(columns, top_names):
                row[column_name] = organism

    with Path(args.output).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["ID", *ordered_columns], delimiter="\t")
        writer.writeheader()
        for sample_id in sorted(sample_rows):
            writer.writerow(sample_rows[sample_id])


if __name__ == "__main__":
    main()
