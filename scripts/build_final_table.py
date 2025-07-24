#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import json
from pathlib import Path
from functools import reduce
from typing import Dict, List, Union

import pandas as pd


# ----------------- кастомные настройки -----------------
RENAME_MAP: Dict[str, str] = {
    # пример: "reads_mapped_percent": "MappedPercent",
}

MODULE_METRICS: Dict[str, List[str]] = {
    "samtools": ["reads_mapped_percent"],
    "bcftools": ["number_of_records", "number_of_SNPs", "number_of_indels"],
    "picard": [
        "MEAN_COVERAGE", "SD_COVERAGE", "MEDIAN_COVERAGE",
        "PCT_1X", "PCT_5X", "PCT_10X", "PCT_30X", "PCT_50X"
    ],
    "fastqc": ["percent_gc", "avg_sequence_length"],
}


# ----------------- argparse helpers -----------------
def csv_list(v: str) -> list[str]:
    return [x.strip() for x in v.split(",") if x.strip()]


def zfill_map(v: str) -> Dict[str, int]:
    out: Dict[str, int] = {}
    if not v:
        return out
    for chunk in v.split(","):
        col, ln = chunk.split(":", 1)
        out[col.strip()] = int(ln)
    return out


# ----------------- форматирование чисел -----------------
def format_with_comma(val: Union[int, float, str, None]) -> str:
    """
    -> '1 234,56' (два знака, запятая вместо точки, '' для NaN/None/пустого).
    Используем не‑разрывный пробел для тысяч, чтобы не конфликтовать с табами/запятой.
    """
    if val == "" or val is None:
        return ""
    try:
        num = float(val)
    except ValueError:
        return str(val)
    return f"{num:,.2f}".replace(",", " ").replace(".", ",")  # U+202F narrow no‑break space


# ----------------- MultiQC → DataFrame -----------------
def multiqc_to_df(path: Union[str, Path]) -> pd.DataFrame:
    p = Path(path)
    opener = gzip.open if p.suffix == ".gz" else open
    with opener(p, "rt", encoding="utf-8") as fh:
        j = json.load(fh)

    data = j.get("report_general_stats_data", {})
    rows: Dict[str, Dict[str, Union[int, float, str]]] = {}

    for module, samples in data.items():
        if module not in MODULE_METRICS:
            continue
        need = MODULE_METRICS[module]
        for sample, metrics in samples.items():
            row = rows.setdefault(sample, {})
            for m in need:
                if m in metrics:
                    col = m if m not in row else f"{module}:{m}"
                    row[col] = metrics[m]

    df = (
        pd.DataFrame.from_dict(rows, orient="index")
        .reset_index()
        .rename(columns={"index": "ID"})
        .astype({"ID": "string"})
        .rename(columns=RENAME_MAP or {})
    )
    return df


# ----------------- Таблицы CSV/TSV -----------------
def read_table(path: Path, str_cols: list[str]) -> pd.DataFrame:
    dtypes = {c: "string" for c in str_cols if c != "ID"}
    df = pd.read_csv(path, sep=None, engine="python", dtype=dtypes, keep_default_na=False)
    first = df.columns[0]
    if first != "ID":
        df = df.rename(columns={first: "ID"})
    df["ID"] = df["ID"].astype("string")
    return df


# ----------------- Merge util -----------------
def merge_tables(dfs: list[pd.DataFrame], how: str = "outer") -> pd.DataFrame:
    return reduce(lambda l, r: pd.merge(l, r, on="ID", how=how), dfs)


# ----------------- main -----------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Собрать таблицу из MultiQC JSON и объединить с другими таблицами."
    )
    p.add_argument("-m", "--multiqc", required=True)
    p.add_argument("-t", "--tables", nargs="*", default=[])
    p.add_argument("-o", "--out", default="FINAL_TABLE.tsv")
    p.add_argument("--join", default="outer", choices=["outer", "inner", "left", "right"])
    p.add_argument("--str-cols", type=csv_list, default=[],
                   help="Колонки читать как строки (сохраняют ведущие нули)")
    p.add_argument("--zfill-cols", type=zfill_map, default={},
                   help="Дополнить слева нулями: col:len,col2:len2")
    p.add_argument("--comma-cols", type=csv_list, default=[],
                   help="Колонки, где точку заменить на запятую и оставить 2 знака")
    return p.parse_args()


def main():
    args = parse_args()

    df_mqc = multiqc_to_df(args.multiqc)
    dfs = [df_mqc]

    # остальные таблицы
    paths: list[Path] = []
    for mask in args.tables:
        paths.extend(sorted(Path().glob(mask)))
    for pth in paths:
        dfs.append(read_table(pth, args.str_cols))

    merged = merge_tables(dfs, how=args.join)

    # zfill
    for col, width in args.zfill_cols.items():
        if col in merged:
            merged[col] = merged[col].astype("string").str.zfill(width)

    # comma‑format
    for col in args.comma_cols:
        if col in merged:
            merged[col] = merged[col].apply(format_with_comma)

    merged.to_csv(args.out, sep="\t", index=False)
    print(f"Готово: {args.out}")


if __name__ == "__main__":
    main()
