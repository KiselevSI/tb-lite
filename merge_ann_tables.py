#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Merge bcftools-query TSVs (POS, REF, ALT, QUAL, ANN) into one table.

Поведение:
- Парсим ANN из КАЖДОЙ таблицы (берём первую запись ANN до запятой, затем делим по '|').
- Для каждого ANN_* собираем значения со всех таблиц для позиции:
    одинаковые -> одно значение
    разные     -> уникальные значения, соединённые ';' в порядке таблиц
- Для каждого входного файла добавляем колонку с ALT, агрегируя повторы ';'.
- Джойним по ключам (по умолчанию: если CHROM есть везде -> CHROM POS REF, иначе POS REF).

Примеры:
    python merge_ann_tables.py -t results2/annotate_vcf/*.tsv -o merged.tsv
    python merge_ann_tables.py -t a.tsv b.tsv -o merged.tsv --keys CHROM POS REF
"""
import argparse
from pathlib import Path
import sys
import pandas as pd
import numpy as np

ANN_FIELDS = [
    "ANN_Allele","ANN_Annotation","ANN_Impact","ANN_Gene_Name","ANN_Gene_ID",
    "ANN_Feature_Type","ANN_Feature_ID","ANN_BioType","ANN_Rank",
    "ANN_HGVS_c","ANN_HGVS_p","ANN_cDNA_pos_len","ANN_CDS_pos_len",
    "ANN_AA_pos_len","ANN_Distance","ANN_Errors",
]

def sample_name_from_path(p: Path) -> str:
    base = p.name
    return base.split('.')[0] if '.' in base else base

def read_table(path: Path) -> pd.DataFrame:
    """
    Ждём колонки POS, REF, ALT, QUAL, ANN (опционально CHROM).
    Если заголовка нет — пытаемся назначить автоматически.
    """
    try:
        df = pd.read_csv(path, sep='\t', dtype=str, na_filter=False)
        df.columns = [str(c).strip() for c in df.columns]
    except Exception as e:
        raise SystemExit(f"[ERROR] cannot read {path}: {e}")

    required = {"POS","REF","ALT","QUAL","ANN"}
    missing = required - set(df.columns)
    if missing:
        # headerless fallback
        df = pd.read_csv(path, sep='\t', dtype=str, na_filter=False, header=None)
        if df.shape[1] >= 6:
            df = df.iloc[:, :6]; df.columns = ["CHROM","POS","REF","ALT","QUAL","ANN"]
        elif df.shape[1] >= 5:
            df = df.iloc[:, :5]; df.columns = ["POS","REF","ALT","QUAL","ANN"]
        else:
            raise SystemExit(f"[ERROR] {path}: expected at least 5 columns (POS REF ALT QUAL ANN)")

    # нормализация
    if "CHROM" in df.columns:
        df["CHROM"] = df["CHROM"].astype(str).str.strip()
    try:
        df["POS"] = pd.to_numeric(df["POS"], errors="raise")
    except Exception:
        pass
    for col in ("ALT","ANN","QUAL","REF"):
        if col in df.columns:
            df[col] = df[col].astype(str).replace({"":".","NA":".","NaN":".","nan":"."})
    return df

def join_unique_order(vals):
    """Соединить уникальные значения, сохраняя порядок появления."""
    seen, out = set(), []
    for v in vals:
        v = str(v)
        if v not in seen:
            seen.add(v); out.append(v)
    return ";".join(out)

def aggregate_alt(df: pd.DataFrame, keys):
    """Группируем по keys и склеиваем ALT ';' с сохранением порядка."""
    return df.groupby(keys, dropna=False, as_index=False)["ALT"].agg(join_unique_order)

def aggregate_ann_rowwise(df: pd.DataFrame, keys) -> pd.DataFrame:
    """
    Для одной таблицы:
      - берём первый непустой ANN и первую QUAL по keys
      - парсим ANN (первая запись до запятой) в ANN_FIELDS
    Возвращает: keys + QUAL + ANN_FIELDS
    """
    def first_nonempty(vals):
        for v in vals:
            if v and v != ".": return v
        return "."
    g = df.groupby(keys, dropna=False, as_index=False).agg(
        ANN=("ANN", first_nonempty),
        QUAL=("QUAL", lambda x: next(iter(x), np.nan))
    )
    ann_series = g["ANN"].fillna(".").replace("", ".").apply(
        lambda s: s.split(",")[0] if s and s!="." else "."
    )
    ann_df = ann_series.str.split("|", n=15, expand=True)
    if ann_df.shape[1] < 16:
        ann_df = ann_df.reindex(columns=range(16), fill_value="")
    ann_df.columns = ANN_FIELDS
    for c in ANN_FIELDS:
        ann_df[c] = ann_df[c].astype(str).str.strip()
    block = pd.concat([g[keys + ["QUAL"]].reset_index(drop=True),
                       ann_df.reset_index(drop=True)], axis=1)
    return block

def aggregate_ann_across_tables(blocks, keys) -> pd.DataFrame:
    """
    На вход: список блоков (keys + QUAL + ANN_FIELDS) по таблицам.
    На выход: keys + QUAL + ANN_FIELDS, где каждый ANN_* — это агрегат по всем таблицам:
      - собираем уникальные значения в порядке таблиц
      - если значение одно, оставляем без ';'
      - QUAL берём первую ненулевую
    """
    frames = []
    for i, b in enumerate(blocks):
        bb = b.copy()
        bb["__src"] = i
        frames.append(bb)
    big = pd.concat(frames, ignore_index=True)
    big = big.sort_values(keys + ["__src"]).reset_index(drop=True)

    def agg_ann(series: pd.Series) -> str:
        vals, seen = [], set()
        for v in series:
            s = str(v)
            if s in ("", "."):  # пропускаем пустое
                continue
            if s not in seen:
                seen.add(s); vals.append(s)
        if not vals:
            return ""
        return vals[0] if len(vals) == 1 else ";".join(vals)

    def agg_qual(series: pd.Series):
        for v in series:
            if pd.notna(v) and str(v) != "":
                return v
        return np.nan

    grouped = big.groupby(keys, dropna=False, as_index=False).agg(
        **({"QUAL": ("QUAL", agg_qual)} | {c: (c, agg_ann) for c in ANN_FIELDS})
    )
    return grouped

def main():
    ap = argparse.ArgumentParser(description="Merge snpEff-ann TSVs; ANN агрегируется по всем таблицам.")
    ap.add_argument("-t","--tables", nargs="+", required=True,
                    help="Input TSVs (уже развёрнутые шеллом, например results2/annotate_vcf/*.tsv)")
    ap.add_argument("-o","--output", required=True, help="Output TSV.")
    ap.add_argument("--keys", nargs="+", default=None,
                    help="Ключи мерджа. Если не указано: CHROM POS REF (если CHROM есть во всех), иначе POS REF.")
    args = ap.parse_args()

    table_paths = [Path(p) for p in args.tables]
    if not table_paths:
        raise SystemExit("[ERROR] -t gave no files (shell glob matched nothing)")

    # читаем все таблицы
    tables = [read_table(p) for p in table_paths]

    # выбираем ключи
    if args.keys is None:
        if all("CHROM" in t.columns for t in tables):
            keys = ["CHROM", "POS", "REF"]
        else:
            keys = ["POS", "REF"]
    else:
        keys = args.keys

    # блоки ANN и ALT по таблицам
    ann_blocks = [aggregate_ann_rowwise(t, keys) for t in tables]
    alt_blocks = [aggregate_alt(t, keys).rename(columns={"ALT": sample_name_from_path(table_paths[i])})
                  for i, t in enumerate(tables)]

    # агрегируем ANN по всем таблицам
    ann_merged = aggregate_ann_across_tables(ann_blocks, keys)

    # приклеиваем sample-колонки ALT
    merged = ann_merged.copy()
    for ab in alt_blocks:
        merged = merged.merge(ab, on=keys, how="outer")

    # финальный порядок колонок
    sample_cols = [sample_name_from_path(p) for p in table_paths]
    for c in ANN_FIELDS:
        if c not in merged.columns:
            merged[c] = ""
    if "QUAL" not in merged.columns:
        merged["QUAL"] = np.nan

    col_order = [c for c in (keys + ["QUAL"] + ANN_FIELDS + sample_cols) if c in merged.columns]
    merged = merged[col_order]

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
