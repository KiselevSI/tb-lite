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


# ----------------- Утилиты для argparse -----------------
def csv_list(value: str) -> List[str]:
    """Парсер 'a,b,c' -> ['a','b','c']"""
    return [x.strip() for x in value.split(",") if x.strip()]


def zfill_map(value: str) -> Dict[str, int]:
    """
    Парсер 'Spol8:8,Barcode:12' -> {'Spol8': 8, 'Barcode': 12}
    """
    out: Dict[str, int] = {}
    if not value:
        return out
    for chunk in value.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if ":" not in chunk:
            raise argparse.ArgumentTypeError(
                f"Формат для --zfill-cols: 'col:len,col2:len2', а не '{chunk}'"
            )
        col, ln = chunk.split(":", 1)
        col = col.strip()
        try:
            ln = int(ln)
        except ValueError:
            raise argparse.ArgumentTypeError(f"Длина должна быть числом: '{chunk}'")
        out[col] = ln
    return out


# ----------------- Аргументы -----------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Собрать таблицу из MultiQC JSON и объединить с другими таблицами по первой колонке."
    )
    p.add_argument("-m", "--multiqc", required=True, help="multiqc_data.json или .json.gz")
    p.add_argument("-t", "--tables", nargs="*", default=[], help="Доп. таблицы (CSV/TSV), можно маски")
    p.add_argument("-o", "--out", default="FINAL_TABLE.tsv", help="Выходной файл (TSV)")
    p.add_argument("--join", default="outer", choices=["outer", "inner", "left", "right"],
                   help="Тип объединения (merge how), по умолчанию outer")
    p.add_argument(
        "--str-cols",
        type=csv_list,
        default=[],
        help="Колонки (через запятую), которые читать как строки (сохранить ведущие нули)."
    )
    p.add_argument(
        "--zfill-cols",
        type=zfill_map,
        default={},
        help="Дополнять слева нулями указанные колонки: col:len,col2:len2"
    )
    return p.parse_args()


# ----------------- Чтение/парсинг MultiQC -----------------
def multiqc_to_df(path: Union[str, Path]) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Файл не найден: {p}")

    opener = gzip.open if p.suffix == ".gz" else open
    with opener(p, "rt", encoding="utf-8") as fh:
        j = json.load(fh)

    data = j.get("report_general_stats_data", {})
    if not isinstance(data, dict):
        raise ValueError("report_general_stats_data имеет неожиданную структуру")

    rows: Dict[str, Dict[str, Union[int, float, str]]] = {}

    for module, samples in data.items():
        if module not in MODULE_METRICS:
            continue
        metrics_needed = MODULE_METRICS[module]
        if not isinstance(samples, dict):
            continue

        for sample, metrics in samples.items():
            if not isinstance(metrics, dict):
                continue
            row = rows.setdefault(sample, {})
            for m in metrics_needed:
                if m in metrics:
                    col = m if m not in row else f"{module}:{m}"
                    # (samtools_3 не берём, поэтому конфликтов тут почти не будет)
                    row[col] = metrics[m]

    df = pd.DataFrame.from_dict(rows, orient="index").reset_index().rename(columns={"index": "ID"})
    # Переименуем по RENAME_MAP
    if RENAME_MAP:
        df = df.rename(columns=RENAME_MAP)
    return df


# ----------------- Чтение прочих таблиц -----------------
def read_table(path: Path, str_cols: List[str]) -> pd.DataFrame:
    """
    Читаем таблицу. Первую колонку переименуем в 'ID'.
    Указанные в str_cols читаем как строки (dtype='string'), остальное autodetect.
    """
    # dtype только для указанных колонок (кроме ID, её приведём вручную)
    dtypes = {c: "string" for c in str_cols if c != "ID"}

    try:
        df = pd.read_csv(
            path,
            sep=None,
            engine="python",
            dtype=dtypes,
            keep_default_na=False
        )
    except Exception:
        df = pd.read_csv(
            path,
            sep="\t",
            dtype=dtypes,
            keep_default_na=False
        )

    # гарантируем, что ID строковый
    first_col = df.columns[0]
    if first_col != "ID":
        df = df.rename(columns={first_col: "ID"})
    df["ID"] = df["ID"].astype("string")
    return df


# ----------------- Merge -----------------
def merge_tables(dfs: List[pd.DataFrame], how: str = "outer") -> pd.DataFrame:
    return reduce(lambda l, r: pd.merge(l, r, on="ID", how=how), dfs)


# ----------------- main -----------------
def main():
    args = parse_args()

    # 1) JSON -> DF
    df_mqc = multiqc_to_df(args.multiqc)
    df_mqc["ID"] = df_mqc["ID"].astype("string")

    # 2) Остальные таблицы
    other_paths: List[Path] = []
    for mask in args.tables:
        other_paths.extend(sorted(Path().glob(mask)))

    dfs = [df_mqc]
    for pth in other_paths:
        dfs.append(read_table(pth, args.str_cols))

    # 3) Merge
    merged = merge_tables(dfs, how=args.join)

    # 4) zfill по запросу
    for col, width in args.zfill_cols.items():
        if col in merged.columns:
            merged[col] = merged[col].astype("string").str.zfill(width)

    # 5) Сохранение
    merged.to_csv(args.out, sep="\t", index=False)
    print(f"Готово: {args.out}")


if __name__ == "__main__":
    main()