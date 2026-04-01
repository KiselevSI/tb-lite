#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd


PICARD_COLUMNS = [
    "MEAN_COVERAGE",
    "SD_COVERAGE",
    "MEDIAN_COVERAGE",
    "PCT_1X",
    "PCT_5X",
    "PCT_10X",
    "PCT_30X",
    "PCT_50X",
]

FINAL_COLUMNS = [
    "ID",
    "number_of_records",
    "number_of_SNPs",
    "number_of_indels",
    "MEAN_COVERAGE",
    "SD_COVERAGE",
    "MEDIAN_COVERAGE",
    "PCT_1X",
    "PCT_5X",
    "PCT_10X",
    "PCT_30X",
    "PCT_50X",
    "reads_mapped_percent",
]


def strip_suffix(name: str, patterns: list[str]) -> str:
    for pattern in patterns:
        stripped = re.sub(pattern, "", name)
        if stripped != name:
            return stripped
    return Path(name).stem


def parse_args():
    parser = argparse.ArgumentParser(
        description="Собрать таблицу основных QC/variant метрик без зависимости от MultiQC."
    )
    parser.add_argument("--wgs", nargs="*", default=[], help="Picard WGS metrics files")
    parser.add_argument("--bcftools", nargs="*", default=[], help="bcftools stats files")
    parser.add_argument("--flagstat", nargs="*", default=[], help="samtools flagstat files")
    parser.add_argument("-o", "--out", default="general.tsv", help="Output TSV path")
    parser.add_argument(
        "--round",
        action="store_true",
        help="Round numeric columns to 2 decimals",
    )
    return parser.parse_args()


def parse_picard(path: Path) -> dict:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    header_idx = next((idx for idx, line in enumerate(lines) if line.startswith("GENOME_TERRITORY\t")), None)
    if header_idx is None or header_idx + 1 >= len(lines):
        return {}

    header = lines[header_idx].split("\t")
    values = lines[header_idx + 1].split("\t")
    row = dict(zip(header, values))
    sample = strip_suffix(
        path.name,
        [
            r"\.CollectWgsMetrics\.coverage_metrics$",
            r"\.wgs_metrics\.txt$",
        ],
    )

    parsed = {"ID": sample}
    for column in PICARD_COLUMNS:
        parsed[column] = row.get(column)
    return parsed


def parse_bcftools(path: Path) -> dict:
    sample = strip_suffix(
        path.name,
        [
            r"\.bcftools_stats\.txt$",
            r"\.bcftools\.txt$",
        ],
    )
    parsed = {"ID": sample}
    mapping = {
        "number of records": "number_of_records",
        "number of SNPs": "number_of_SNPs",
        "number of indels": "number_of_indels",
    }

    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("SN\t"):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        key = parts[2].rstrip(":")
        if key in mapping:
            parsed[mapping[key]] = parts[3]
    return parsed


def parse_flagstat(path: Path) -> dict:
    sample = strip_suffix(
        path.name,
        [
            r"\.flagstat$",
            r"\.flagstat\.txt$",
        ],
    )
    parsed = {"ID": sample}
    match = re.search(r"mapped \(([0-9.]+)%", path.read_text(encoding="utf-8"))
    if match:
        parsed["reads_mapped_percent"] = match.group(1)
    return parsed


def main():
    args = parse_args()

    rows: dict[str, dict] = {}

    for entry in args.wgs:
        parsed = parse_picard(Path(entry))
        if parsed:
            rows.setdefault(parsed["ID"], {}).update(parsed)

    for entry in args.bcftools:
        parsed = parse_bcftools(Path(entry))
        if parsed:
            rows.setdefault(parsed["ID"], {}).update(parsed)

    for entry in args.flagstat:
        parsed = parse_flagstat(Path(entry))
        if parsed:
            rows.setdefault(parsed["ID"], {}).update(parsed)

    df = pd.DataFrame(rows.values()) if rows else pd.DataFrame(columns=["ID"])

    for column in FINAL_COLUMNS:
        if column not in df.columns:
            df[column] = pd.NA

    df = df[FINAL_COLUMNS].sort_values("ID", kind="stable")

    if args.round:
        for column in FINAL_COLUMNS:
            if column == "ID":
                continue
            df[column] = pd.to_numeric(df[column], errors="coerce").round(2)

    df.to_csv(
        args.out,
        sep="\t",
        index=False,
        float_format="%.2f" if args.round else None,
        decimal=",",
    )


if __name__ == "__main__":
    main()
