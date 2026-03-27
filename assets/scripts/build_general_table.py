#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import json
from pathlib import Path
from typing import Dict, List, Union

import pandas as pd


MODULE_METRICS: Dict[str, List[str]] = {
    "samtools": ["reads_mapped_percent"],
    "bcftools": ["number_of_records", "number_of_SNPs"],
    "picard": [
        "MEAN_COVERAGE",
        "MEDIAN_COVERAGE",
        "PCT_1X",
        "PCT_5X",
        "PCT_10X",
        "PCT_30X",
        "PCT_50X",
    ],
    "fastqc": ["percent_gc", "avg_sequence_length"],
}

FINAL_COLUMNS = [
    "ID",
    "number_of_records",
    "number_of_SNPs",
    "MEAN_COVERAGE",
    "MEDIAN_COVERAGE",
    "PCT_1X",
    "PCT_5X",
    "PCT_10X",
    "PCT_30X",
    "PCT_50X",
    "reads_mapped_percent",
    "percent_gc",
    "avg_sequence_length",
]


def multiqc_to_df(path: Union[str, Path]) -> pd.DataFrame:
    p = Path(path)
    opener = gzip.open if p.suffix == ".gz" else open

    with opener(p, "rt", encoding="utf-8") as fh:
        data = json.load(fh)

    rows: Dict[str, Dict[str, Union[int, float, str]]] = {}

    for module_name, samples in data.get("report_general_stats_data", {}).items():
        if module_name not in MODULE_METRICS:
            continue

        for sample_name, metrics in samples.items():
            row = rows.setdefault(sample_name, {})
            for metric in MODULE_METRICS[module_name]:
                if metric in metrics:
                    row[metric] = metrics[metric]

    df = (
        pd.DataFrame.from_dict(rows, orient="index")
        .reset_index()
        .rename(columns={"index": "ID"})
    )

    df["ID"] = df["ID"].astype("string")

    for col in FINAL_COLUMNS:
        if col not in df.columns:
            df[col] = pd.NA

    return df[FINAL_COLUMNS]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Извлечь выбранные метрики из MultiQC JSON в TSV."
    )
    parser.add_argument(
        "-m", "--multiqc",
        required=True,
        help="Путь к multiqc_data.json или multiqc_data.json.gz"
    )
    parser.add_argument(
        "-o", "--out",
        default="FINAL_TABLE.tsv",
        help="Выходной TSV-файл"
    )
    parser.add_argument(
        "--round",
        action="store_true",
        help="Округлить числовые колонки до 2 знаков"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    df = multiqc_to_df(args.multiqc)

    if args.round:
        for col in FINAL_COLUMNS:
            if col != "ID":
                df[col] = pd.to_numeric(df[col], errors="coerce").round(2)

    df.to_csv(
        args.out,
        sep="\t",
        index=False,
        float_format="%.2f" if args.round else None,
        decimal=","
    )

    print(f"Готово: {args.out}")


if __name__ == "__main__":
    main()