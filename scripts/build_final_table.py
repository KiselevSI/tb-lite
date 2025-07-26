#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import json
import re
from pathlib import Path
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
        "PCT_1X", "PCT_5X", "PCT_10X", "PCT_30X", "PCT_50X",
    ],
    "fastqc": ["percent_gc", "avg_sequence_length"],
}


# ----------------- argparse helpers -----------------
def csv_list(v: str) -> list[str]:
    return [x.strip() for x in v.split(",") if x.strip()]


def zfill_map(v: str) -> Dict[str, int]:
    out: Dict[str, int] = {}
    if v:
        for chunk in v.split(","):
            col, ln = chunk.split(":", 1)
            out[col.strip()] = int(ln)
    return out


def _safe_name(*parts: Union[str, None]) -> str:
    s = "_".join([p for p in parts if isinstance(p, str) and p and p.strip()])
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^0-9A-Za-z_]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "col"


# ----------------- MultiQC → DataFrame -----------------
def multiqc_to_df(path: Union[str, Path]) -> pd.DataFrame:
    p = Path(path)
    opener = gzip.open if p.suffix == ".gz" else open
    with opener(p, "rt", encoding="utf-8") as fh:
        j = json.load(fh)

    rows: Dict[str, Dict[str, Union[int, float, str]]] = {}
    for mod, samples in j.get("report_general_stats_data", {}).items():
        if mod not in MODULE_METRICS:
            continue
        for sample, metrics in samples.items():
            row = rows.setdefault(sample, {})
            for m in MODULE_METRICS[mod]:
                if m in metrics:
                    col = m if m not in row else f"{mod}:{m}"
                    row[col] = metrics[m]

    return (
        pd.DataFrame.from_dict(rows, orient="index")
        .reset_index()
        .rename(columns={"index": "ID"})
        .astype({"ID": "string"})
        .rename(columns=RENAME_MAP or {})
    )


# ----------------- DR Excel → плоская таблица + карты уровней -----------------
def dr_to_flat_with_maps(path: Union[str, Path]) -> tuple[pd.DataFrame, Dict[str, str], Dict[str, str]]:
    """
    Вернуть (dr_flat, top_map, sub_map).

    - Читаем Excel с двухуровневой шапкой: header=[0,1].
    - ID берём из колонки 'Sample' (может быть на любом уровне).
    - Сплющиваем имена в '<Top>_<Sub>'; для 'Other Variants' — фиксированное 'Other_Variants'.
    - top_map: {flat -> верхний уровень (препарат / Unnamed...)}
    - sub_map: {flat -> нижний уровень (gene_name / Mutation / Freq / Other Variants ...)}
    """
    raw = pd.read_excel(path, header=[0, 1])  # MultiIndex columns

    # Ищем колонку Sample
    sample_cols = [
        c for c in raw.columns
        if isinstance(c, tuple)
        and any(isinstance(c[i], str) and c[i].strip().lower() == "sample" for i in range(len(c)))
    ]
    if not sample_cols:
        raise ValueError("В DR-файле не найдена колонка 'Sample'.")
    sample_col = sample_cols[0]

    id_series = raw[sample_col].astype("string").str.strip()

    # ВАЖНО: не удаляем all-NaN колонки — они нужны для сохранения структуры!
    # raw2 = raw.dropna(axis=1, how="all")  # НЕЛЬЗЯ

    cols: list[str] = ["ID"]
    top_map: Dict[str, str] = {}
    sub_map: Dict[str, str] = {}
    series_list: list[pd.Series] = [id_series]

    for col in raw.columns:
        if col == sample_col:
            continue

        if isinstance(col, tuple):
            top = col[0] if len(col) > 0 else ""
            sub = col[1] if len(col) > 1 else ""
        else:
            top, sub = str(col), ""

        top_s = "" if top is None else str(top)
        sub_s = "" if sub is None else str(sub)
        sub_norm = sub_s.strip().lower().replace(" ", "_")

        # Спец-имя для Other Variants
        if sub_norm == "other_variants":
            flat = "Other_Variants"
        else:
            flat = _safe_name(top_s or None, sub_s or None) or "col"

        # уникализируем
        base = flat
        k = 1
        while flat in cols:
            k += 1
            flat = f"{base}__{k}"

        cols.append(flat)
        top_map[flat] = top_s
        sub_map[flat] = sub_s
        series_list.append(raw[col])

    dr_flat = pd.concat(series_list, axis=1)
    dr_flat.columns = cols
    dr_flat["ID"] = dr_flat["ID"].astype("string")
    dr_flat = dr_flat[dr_flat["ID"].notna() & (dr_flat["ID"] != "")]
    return dr_flat.reset_index(drop=True), top_map, sub_map


# ----------------- Таблицы CSV/TSV -----------------
def read_table(path: Path, str_cols: list[str]) -> pd.DataFrame:
    dtypes = {c: "string" for c in str_cols if c != "ID"}
    df = pd.read_csv(
        path,
        sep=None,            # autodetect
        engine="python",     # требуется при sep=None
        dtype=dtypes,
        keep_default_na=False,
    )
    if df.columns[0] != "ID":
        df = df.rename(columns={df.columns[0]: "ID"})
    df["ID"] = df["ID"].astype("string")
    return df


# ----------------- Merge util -----------------
def merge_tables(dfs: list[pd.DataFrame], how: str = "outer") -> pd.DataFrame:
    """
    Объединяем список DataFrame по ID.
    Переводим в индекс и последовательно конкатенируем столбцы:
    pd.concat([base, r], axis=1, join=how).

    Для одноуровневых пересечений имён добавляем *_x / *_y.
    """
    if not dfs:
        return pd.DataFrame(columns=["ID"])

    frames = [df.set_index("ID") for df in dfs]
    base = frames[0].copy()
    suffixed_once: set[str] = set()

    for r in frames[1:]:
        base_one = {c for c in base.columns if not isinstance(c, tuple)}
        r_one = {c for c in r.columns if not isinstance(c, tuple)}
        overlap = sorted(base_one.intersection(r_one))

        if overlap:
            # base: *_x
            rename_base: Dict[str, str] = {}
            for c in overlap:
                if c not in suffixed_once and c in base.columns:
                    new = f"{c}_x"
                    k = 1
                    while new in base.columns:
                        k += 1
                        new = f"{c}_x{k}"
                    rename_base[c] = new
                    suffixed_once.add(c)
            if rename_base:
                base = base.rename(columns=rename_base)

            # r: *_y
            rename_r: Dict[str, str] = {}
            for c in overlap:
                if c in r.columns:
                    new = f"{c}_y"
                    k = 1
                    while new in r.columns or new in base.columns:
                        k += 1
                        new = f"{c}_y{k}"
                    rename_r[c] = new
            if rename_r:
                r = r.rename(columns=rename_r)

        base = pd.concat([base, r], axis=1, join=how)

    return base.reset_index()


# ----------------- main -----------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Собрать итоговую таблицу из MultiQC JSON и других файлов."
    )
    p.add_argument("-m", "--multiqc", required=True,
                   help="multiqc_data.json или .json.gz")
    p.add_argument("-t", "--tables", nargs="*", default=[],
                   help="Дополнительные таблицы (маски путей допускаются)")
    p.add_argument("--dr", default=None,
                   help="Excel-файл с DR. Сливается по колонке Sample → ID")
    p.add_argument("-o", "--out", default="FINAL_TABLE.xlsx",
                   help="Имя выходного файла (.xlsx или .tsv)")
    p.add_argument("--join", choices=["outer", "inner", "left", "right"],
                   default="outer", help="Тип объединения")
    p.add_argument("--str-cols", type=csv_list, default=[],
                   help="Колонки читать как строки (сохраняют ведущие нули)")
    p.add_argument("--zfill-cols", type=zfill_map, default={},
                   help="Дополнить слева нулями: col:len,col2:len2")
    p.add_argument("--round-cols", type=csv_list, default=[],
                   help="Округлить до 2 знаков (при этом остаются числа)")
    return p.parse_args()


def main():
    a = parse_args()

    # 1) База без DR
    base_dfs: list[pd.DataFrame] = [multiqc_to_df(a.multiqc)]
    for mask in a.tables:
        for p in sorted(Path().glob(mask)):
            base_dfs.append(read_table(p, a.str_cols))
    base = merge_tables(base_dfs, how=a.join)

    # 2) DR
    top_map: Dict[str, str] = {}
    sub_map: Dict[str, str] = {}
    if a.dr:
        dr_flat, top_map, sub_map = dr_to_flat_with_maps(a.dr)
        base = pd.merge(base, dr_flat, on="ID", how=a.join)

    merged = base

    # 3) Постобработка значений
    for col, width in a.zfill_cols.items():
        if col in merged:
            merged[col] = merged[col].astype("string").str.zfill(width)
    for col in a.round_cols:
        if col in merged:
            merged[col] = pd.to_numeric(merged[col], errors="coerce").round(2)

    def try_numeric(val):
        try:
            return pd.to_numeric(val)
        except (ValueError, TypeError):
            return val

    for col in merged.columns:
        if col not in a.str_cols and merged[col].dtype == "object":
            merged[col] = merged[col].apply(try_numeric)

    # 4) Вывод
    write_df = merged.copy()
    has_dr = bool(top_map)

    if a.out.lower().endswith(".xlsx"):
        with pd.ExcelWriter(a.out, engine="xlsxwriter") as writer:
            if has_dr:
                # Строка 0 (верхняя): препарат только над *_Mutation
                # Строка 1: подзаголовки
                top_row: list[str] = []
                sub_row: list[str] = []

                for col in write_df.columns:
                    top = top_map.get(col, "")
                    sub = sub_map.get(col, "")
                    sub_norm = re.sub(r"\s+", "_", str(sub).strip().lower())

                    # верхняя строка
                    if top and sub_norm == "mutation":
                        top_row.append(top)
                    else:
                        top_row.append("")

                    # строка подзаголовков
                    if col in sub_map:
                        if sub_norm == "gene_name":
                            sub_row.append("gene_name")
                        elif sub_norm == "mutation":
                            sub_row.append("Mutation")
                        elif sub_norm == "freq":
                            sub_row.append("Freq")
                        elif sub_norm == "other_variants":
                            sub_row.append("Other_Variants")
                        else:
                            sub_row.append(sub)
                    else:
                        sub_row.append(col)

                # пишем две строки шапки вручную
                write_df.to_excel(writer, sheet_name="Data", startrow=2, index=False, header=False)
                ws = writer.sheets["Data"]
                ws.write_row(0, 0, top_row)
                ws.write_row(1, 0, sub_row)
            else:
                write_df.to_excel(writer, sheet_name="Data", index=False)

            # числовой формат
            try:
                num_fmt = writer.book.add_format({"num_format": "#,##0.00"})
                ws = writer.sheets["Data"]
                for idx, col in enumerate(write_df.columns):
                    if pd.api.types.is_numeric_dtype(write_df[col]):
                        ws.set_column(idx, idx, None, num_fmt)
            except Exception:
                pass
    else:
        # TSV
        if has_dr:
            top_row: list[str] = []
            sub_row: list[str] = []
            for col in write_df.columns:
                top = top_map.get(col, "")
                sub = sub_map.get(col, "")
                sub_norm = re.sub(r"\s+", "_", str(sub).strip().lower())

                if top and sub_norm == "mutation":
                    top_row.append(top)
                else:
                    top_row.append("")

                if col in sub_map:
                    if sub_norm == "gene_name":
                        sub_row.append("gene_name")
                    elif sub_norm == "mutation":
                        sub_row.append("Mutation")
                    elif sub_norm == "freq":
                        sub_row.append("Freq")
                    elif sub_norm == "other_variants":
                        sub_row.append("Other_Variants")
                    else:
                        sub_row.append(sub)
                else:
                    sub_row.append(col)

            with open(a.out, "w", encoding="utf-8") as fh:
                fh.write("\t".join(top_row) + "\n")
                fh.write("\t".join(sub_row) + "\n")
                write_df.to_csv(
                    fh,
                    sep="\t",
                    index=False,
                    header=False,   # шапку записали вручную
                    float_format="%.2f",
                    decimal=",",
                )
        else:
            write_df.to_csv(
                a.out,
                sep="\t",
                index=False,
                float_format="%.2f",
                decimal=",",
            )

    print(f"Готово: {a.out}")


if __name__ == "__main__":
    main()
