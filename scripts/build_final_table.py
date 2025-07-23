#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import json
from pathlib import Path
from functools import reduce
from typing import Dict, List, Union

import pandas as pd


# ======= НАСТРОЙКИ ПОЛЕЙ (впиши свои имена) =======
RENAME_MAP: Dict[str, str] = {
    "number_of_records": "Vars",
    "number_of_SNPs": "SNP",
    "number_of_indels": "Indel",
    # Пример:
    # "reads_mapped_percent": "Mapped_%"
}


# Какие метрики берём из каких модулей
MODULE_METRICS: Dict[str, List[str]] = {
    "samtools":   ["reads_mapped_percent"],
    # "samtools_3": ["reads_mapped_percent"],  # <-- НЕ ИСПОЛЬЗУЕМ
    "bcftools":   ["number_of_records", "number_of_SNPs", "number_of_indels"],
    "picard":     ["MEAN_COVERAGE", "MEDIAN_COVERAGE", "SD_COVERAGE", 
                   "PCT_1X", "PCT_5X", "PCT_10X", "PCT_30X", "PCT_50X"],
    "fastqc":     ["percent_gc", "avg_sequence_length", "percent_duplicates"],
}


# ----------------- Аргументы -----------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Собрать таблицу из MultiQC JSON и объединить с другими таблицами по первой колонке."
    )
    p.add_argument("-m", "--multiqc", required=True, help="multiqc_data.json или .json.gz")
    p.add_argument("-t", "--tables", nargs="*", default=[], help="Доп. таблицы (CSV/TSV), можно маски")
    p.add_argument("-o", "--out", default="FINAL.tsv", help="Выходной TSV")
    p.add_argument("--join", default="outer", choices=["outer", "inner", "left", "right"],
                   help="Тип merge (по умолчанию outer)")
    return p.parse_args()


# ----------------- Чтение/парсинг MultiQC -----------------
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
        needed = MODULE_METRICS[module]
        if not isinstance(samples, dict):
            continue

        for sample, metrics in samples.items():
            if not isinstance(metrics, dict):
                continue
            row = rows.setdefault(sample, {})
            for m in needed:
                if m in metrics:
                    # если метрика уже была, добавим префикс модуля, но samtools_3 мы не берём
                    col = m if col_free(row, m) else f"{module}:{m}"
                    row[col] = metrics[m]

    df = pd.DataFrame.from_dict(rows, orient="index").reset_index().rename(columns={"index": "ID"})
    return df


def col_free(row: Dict[str, Union[int, float, str]], name: str) -> bool:
    """True если такого поля ещё нет в строке."""
    return name not in row


# ----------------- Чтение прочих таблиц -----------------
def read_table(path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep=None, engine="python")
    except Exception:
        df = pd.read_csv(path, sep="\t")
    first_col = df.columns[0]
    if first_col != "ID":
        df = df.rename(columns={first_col: "ID"})
    return df


# ----------------- Merge -----------------
def merge_tables(dfs: List[pd.DataFrame], how: str = "outer") -> pd.DataFrame:
    return reduce(lambda l, r: pd.merge(l, r, on="ID", how=how), dfs)


# ----------------- main -----------------
def main():
    args = parse_args()

    # 1) JSON -> DF
    df_mqc = multiqc_to_df(args.multiqc)

    # 2) Переименования (заполни RENAME_MAP выше)
    if RENAME_MAP:
        df_mqc = df_mqc.rename(columns=RENAME_MAP)

    # 3) Остальные таблицы
    other_paths: List[Path] = []
    for mask in args.tables:
        other_paths.extend(sorted(Path().glob(mask)))

    dfs = [df_mqc]
    for p in other_paths:
        dfs.append(read_table(p))

    # 4) Merge
    merged = merge_tables(dfs, how=args.join)

    # 5) Save
    merged.to_csv(args.out, sep=",", index=False)
    print(f"Готово: {args.out}")


if __name__ == "__main__":
    main()
