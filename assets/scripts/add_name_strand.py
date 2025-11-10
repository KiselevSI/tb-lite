#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import re
import pandas as pd

ID_COL = "ANN[*].GENEID"

def read_table(path: str) -> pd.DataFrame:
    """Read TSV/TSV.GZ with UTF-8, everything as str."""
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            dtype=str,
            na_filter=False,
            compression="infer",
            engine="python"
        )
        return df
    except Exception as e:
        print(f"[ERROR] Failed to read {path}: {e}", file=sys.stderr)
        sys.exit(1)

def split_ids(val: str):
    """Split ANN[*].GENEID into tokens by -, ;, , or | . Keeps order, strips spaces."""
    if val is None:
        return []
    s = str(val).strip()
    if not s:
        return []
    # Разделители: дефис, точка с запятой, запятая, вертикальная черта
    tokens = re.split(r"\s*[-;,|]\s*", s)
    return [t for t in tokens if t]

def main():
    ap = argparse.ArgumentParser(
        description="Добавить 'name' и 'strand' из NCBI feature_table в таблицу по locus_tag/ANN[*].GENEID. Поддерживает несколько ID в одной ячейке."
    )
    ap.add_argument("-i", "--input", required=True, help=f"Входной TSV с колонкой {ID_COL}")
    ap.add_argument("-t", "--feature-table", required=True, help="NCBI feature_table.txt(.gz)")
    ap.add_argument("-o", "--output", required=True, help="Путь для выходного TSV")
    args = ap.parse_args()

    # Read user table
    user_df = read_table(args.input)
    if ID_COL not in user_df.columns:
        print(f"[ERROR] Во входной таблице нет колонки '{ID_COL}'.", file=sys.stderr)
        sys.exit(2)

    # Read feature_table and filter to CDS
    ft = read_table(args.feature_table)
    required_cols = {"# feature", "locus_tag", "name", "strand"}
    missing = [c for c in required_cols if c not in ft.columns]
    if missing:
        print(f"[ERROR] В feature_table отсутствуют колонки: {', '.join(missing)}", file=sys.stderr)
        sys.exit(3)

    cds = (
        ft[ft["# feature"] == "CDS"]
        .loc[:, ["locus_tag", "name", "strand"]]
        .dropna(subset=["locus_tag"])
        .drop_duplicates(subset=["locus_tag"], keep="first")
    )

    # Build maps locus_tag -> name/strand
    name_map = dict(zip(cds["locus_tag"], cds["name"]))
    strand_map = dict(zip(cds["locus_tag"], cds["strand"]))

    # Compute name/strand per row, supporting multiple IDs in ANN[*].GENEID
    ids_series = user_df[ID_COL].map(split_ids)

    def join_mapped(tokens, mapping):
        # Для каждого токена берём значение из mapping (или ''), затем объединяем через ';'
        return ";".join([mapping.get(tok, "") for tok in tokens]) if tokens else ""

    user_df["name"] = ids_series.apply(lambda toks: join_mapped(toks, name_map))
    user_df["strand"] = ids_series.apply(lambda toks: join_mapped(toks, strand_map))

    # Вставим колонки сразу после ANN[*].GENEID
    cols = list(user_df.columns)
    # Переместим name/strand, если уже существуют
    for c in ("name", "strand"):
        if c in cols:
            cols.remove(c)
    try:
        idx = cols.index(ID_COL)
    except ValueError:
        print(f"[ERROR] Не удалось найти колонку '{ID_COL}' после обработки.", file=sys.stderr)
        sys.exit(4)

    new_cols = cols[: idx + 1] + ["name", "strand"] + cols[idx + 1 :]
    user_df = user_df[new_cols]

    # Статистика соответствий: строк, где найден хотя бы один токен в feature_table
    matched_rows = ids_series.apply(lambda toks: any(tok in name_map for tok in toks)).sum()
    total_rows = len(user_df)

    # Write output
    try:
        user_df.to_csv(args.output, sep="\t", index=False)
    except Exception as e:
        print(f"[ERROR] Не удалось записать файл {args.output}: {e}", file=sys.stderr)
        sys.exit(5)

    print(
        f"[OK] Записано: {args.output}\n"
        f"Всего строк: {total_rows}\n"
        f"Строк с найденными соответствиями (CDS): {matched_rows}",
        file=sys.stderr,
    )

if __name__ == "__main__":
    pd.options.mode.chained_assignment = None
    main()
