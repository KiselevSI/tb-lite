#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import re
import pandas as pd
from typing import List, Optional

DEFAULT_COLUMNS = [
    "ANN[*].EFFECT",
    "ANN[*].IMPACT",
    "ANN[*].GENE",
    "ANN[*].GENEID",
    "ANN[*].FEATUREID",
    "ANN[*].BIOTYPE",
    "ANN[*].RANK",
    "ANN[*].HGVS_C",
    "ANN[*].HGVS_P",
    "ANN[*].CDS_POS",
    "ANN[*].AA_POS",
    "name",
    "strand",
]

ONLY_SEPS_RE = re.compile(r"^\s*[;,]+\s*$")

def read_table(path: str) -> pd.DataFrame:
    try:
        return pd.read_csv(
            path,
            sep="\t",
            dtype=str,
            na_filter=False,
            compression="infer",
            engine="python",
        )
    except Exception as e:
        print(f"[ERROR] Failed to read {path}: {e}", file=sys.stderr)
        sys.exit(1)

def write_table(df: pd.DataFrame, path: str):
    try:
        df.to_csv(path, sep="\t", index=False)
    except Exception as e:
        print(f"[ERROR] Failed to write {path}: {e}", file=sys.stderr)
        sys.exit(2)

def auto_sep(s: str, forced: Optional[str]) -> str:
    if forced:
        return forced
    has_sc = ";" in s
    has_cm = "," in s
    if has_sc and not has_cm:
        return ";"
    if has_cm and not has_sc:
        return ","
    return ","

def dedup_cell(value: str, forced_sep: Optional[str] = None) -> str:
    """
    Правила:
      - Пустые ячейки (None или "") не трогаем.
      - Если ячейка содержит только разделители ';' и/или ',' (и пробелы) — делаем пустой.
      - Иначе: сплит по ';' или ',', убрать дубли с сохранением порядка, склеить выбранным разделителем.
    """
    # 1) пустые оставляем как есть
    if value is None or value == "":
        return value

    s = str(value)

    # 2) только разделители -> делаем пусто
    if ONLY_SEPS_RE.match(s):
        return ""

    # 3) обычная обработка
    tokens = re.split(r"\s*[;,]\s*", s)
    seen = set()
    uniq = []
    for tok in tokens:
        tok = tok.strip()
        if not tok:
            continue
        if tok not in seen:
            seen.add(tok)
            uniq.append(tok)

    # если после очистки ничего не осталось — делаем пусто
    if not uniq:
        return ""

    sep = auto_sep(s, forced_sep)
    return sep.join(uniq)

def process(df: pd.DataFrame, columns: List[str], forced_sep: Optional[str] = None) -> pd.DataFrame:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        print(f"[WARN] Пропущены несуществующие колонки: {', '.join(missing)}", file=sys.stderr)
    for col in columns:
        if col in df.columns:
            df[col] = df[col].map(lambda v: dedup_cell(v, forced_sep))
    return df

def main():
    ap = argparse.ArgumentParser(
        description="Убрать дубликаты внутри ячеек (разделители ',' и ';'). Пустые ячейки не трогаются. Ячейки с одними разделителями становятся пустыми."
    )
    ap.add_argument("-i", "--input", required=True, help="Входной TSV/TSV.GZ")
    ap.add_argument("-o", "--output", required=True, help="Выходной TSV")
    ap.add_argument(
        "-c", "--columns",
        nargs="+",
        default=DEFAULT_COLUMNS,
        help="Список колонок для обработки (по умолчанию — ANN[*] + name,strand)."
    )
    ap.add_argument(
        "--sep",
        choices=[",", ";"],
        default=None,
        help="Принудительный разделитель результата (по умолчанию авто-выбор по содержимому ячейки)."
    )
    args = ap.parse_args()

    df = read_table(args.input)
    df = process(df, args.columns, args.sep)
    write_table(df, args.output)

    present = [c for c in args.columns if c in df.columns]
    print(
        f"[OK] Обработано колонок: {len(present)} / {len(args.columns)}. "
        f"Вывод записан: {args.output}",
        file=sys.stderr
    )

if __name__ == "__main__":
    pd.options.mode.chained_assignment = None
    main()
