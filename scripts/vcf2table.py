#!/usr/bin/env python3
import argparse
import os
import re
import sys
import pandas as pd
import pysam

# ---------------------------
# Работа с feature_table.tsv
# ---------------------------

def load_feature_tables(file_path: str):
    df = pd.read_csv(file_path, sep="\t")
    feature_table = df[["symbol", "locus_tag"]].drop_duplicates().dropna()
    name_table = df[["name", "locus_tag"]].drop_duplicates()
    return feature_table, name_table

def get_locus_tag(symbol_value: str, feature_table: pd.DataFrame) -> str:
    if not symbol_value:
        return "None"
    hits = feature_table[feature_table["symbol"] == symbol_value]["locus_tag"].values
    return hits[0] if hits.size > 0 else symbol_value

def get_name_tag(locus_tag_value: str, name_table: pd.DataFrame) -> str:
    if not locus_tag_value:
        return "None"
    parts = locus_tag_value.split("-")
    out = []
    for p in parts:
        if not p:
            out.append("None")
            continue
        names = name_table[name_table["locus_tag"] == p]["name"].values
        if len(names) == 0:
            out.append("None")
        elif len(names) == 1:
            out.append(str(names[0]) if names[0] is not None else "None")
        else:
            out.append(str(names[1]) if names[1] is not None else "None")
    return ",".join(out)

def rename_gene_id(gene_name: str, gene_id: str, feature_table: pd.DataFrame) -> str:
    if not gene_id:
        return "None"
    gene_id_parts = gene_id.split("-")
    gene_name_parts = (gene_name or "").split("-")
    while len(gene_name_parts) < len(gene_id_parts):
        gene_name_parts.append("")
    for i, _ in enumerate(gene_id_parts):
        if i < len(gene_name_parts):
            gene_name_parts[i] = get_locus_tag(gene_name_parts[i], feature_table)
    return "-".join(gene_name_parts)

# ---------------------------
# Чтение VCF
# ---------------------------

def read_vcf_pos_ref_alt(path: str) -> pd.DataFrame:
    """Считать только POS, REF, ALT (как 'Allele' = объединение ALT через запятую)."""
    rows = []
    with pysam.VariantFile(path) as vcf:
        for rec in vcf.fetch():
            pos = rec.pos
            ref = rec.ref
            alts = ",".join(rec.alts) if rec.alts else ""
            rows.append((pos, ref, alts))
    return pd.DataFrame(rows, columns=["POS", "REF", "Allele"])

def parse_ann_list(info_ann):
    """Вернуть список строк аннотаций (каждая — '|'‑разделённая)."""
    if info_ann is None:
        return []
    # pysam возвращает tuple/列表; приводим к списку строк
    if isinstance(info_ann, (list, tuple)):
        return [str(x) for x in info_ann]
    return [str(info_ann)]

def process_annotated_record(rec, feature_table, name_table):
    """
    Разбор первой аннотации ANN, при множественных — собираем сводные значения через запятую.
    Формат ANN описан в документации SnpEff.
    """
    anns = parse_ann_list(rec.info.get("ANN"))
    # 16 полей по стандарту; заполним пустыми
    first = [""] * 16
    if anns:
        first_split = anns[0].split("|")
        for i in range(min(16, len(first_split))):
            first[i] = first_split[i]

    allele = ",".join(rec.alts) if rec.alts else ""

    if len(rec.alts or []) <= 1:
        temp_annotation       = first[1]
        temp_put_impact       = first[2]
        temp_feature_type     = first[5]
        temp_hgvsc            = first[9]
        temp_hgvsp            = first[10]
    else:
        ann_list, imp_list, ftype_list, hgvsc_list, hgvsp_list = [], [], [], [], []
        for a in anns:
            s = a.split("|")
            ann_list.append(s[1])
            imp_list.append(s[2])
            ftype_list.append(s[5])
            hgvsc_list.append(s[9])
            hgvsp_list.append(s[10])
        temp_annotation   = ",".join(ann_list)
        temp_put_impact   = ",".join(imp_list)
        temp_feature_type = ",".join(ftype_list)
        temp_hgvsc        = ",".join(hgvsc_list)
        temp_hgvsp        = ",".join(hgvsp_list)

    gene_name = first[3] if len(first) > 3 else ""
    gene_id_raw = first[4] if len(first) > 4 else ""

    try:
        gene_id = rename_gene_id(gene_name, gene_id_raw, feature_table)
    except Exception:
        gene_id = gene_id_raw or "None"

    try:
        name_tag = get_name_tag(gene_id, name_table)
    except Exception:
        name_tag = "None"

    row = [
        rec.pos,
        rec.ref,
        allele,
        temp_annotation,
        temp_put_impact,
        gene_name,
        gene_id,
        name_tag,
        temp_feature_type,
        first[7],
        temp_hgvsc,
        temp_hgvsp,
        first[11],
        first[12],
        first[13],
        first[15],
    ]
    return row

def extract_annotations(annotated_vcf: str, feature_table_path: str) -> pd.DataFrame:
    feature_table, name_table = load_feature_tables(feature_table_path)

    columns = [
        "POS", "REF", "Allele", "Annotation", "Putative_impact", "Gene Name", "Gene ID", "name",
        "Feature type", "Transcript biotype", "HGVS.c", "HGVS.p", "cDNA_position / cDNA_len",
        "CDS_position / CDS_len", "Protein_position / Protein_len", "Errors"
    ]
    rows = []
    with pysam.VariantFile(annotated_vcf) as vcf:
        for rec in vcf.fetch():
            try:
                rows.append(process_annotated_record(rec, feature_table, name_table))
            except Exception as e:
                # Минимальная защитная запись
                rows.append([rec.pos, rec.ref, "", "", "", "", "", "None", "", "", "", "", "", "", "", f"Error: {e}"])
    df = pd.DataFrame(rows, columns=columns).sort_values(by="POS").reset_index(drop=True)
    return df

# ---------------------------
# Добавление колонок по VCF‑списку
# ---------------------------

def update_with_vcfs(df: pd.DataFrame, vcfs: list) -> pd.DataFrame:
    df = df.copy()
    df.set_index("POS", inplace=True)

    for vcf_path in vcfs:
        try:
            short = os.path.basename(vcf_path)
            short = re.sub(r"\.(vcf|vcf\.gz|bcf)(\.tbi)?$", "", short, flags=re.IGNORECASE)
            s = read_vcf_pos_ref_alt(vcf_path).set_index("POS")["Allele"]
            df[short] = s
        except Exception as e:
            sys.stderr.write(f"Предупреждение: не удалось обработать {vcf_path}: {e}\n")

    df.reset_index(inplace=True)
    return df

# ---------------------------
# CLI
# ---------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Извлечь аннотации из готового аннотированного VCF и добавить аллели из списка VCF."
    )
    parser.add_argument("-a", "--annotated", required=True, help="Итоговый аннотированный VCF (с полем ANN).")
    parser.add_argument("-v", "--vcfs", nargs="+", required=True, help="Список VCF файлов для добавления столбцов аллелей.")
    parser.add_argument("-t", "--table", required=True, help="Путь к feature_table.tsv (symbol, locus_tag, name).")
    parser.add_argument("-o", "--output", required=True, help="Путь к выходному .xlsx файлу.")
    args = parser.parse_args()

    # Шаг 1: извлечь аннотации из аннотированного VCF
    df = extract_annotations(args.annotated, args.table)

    # Шаг 2: добавить столбцы по каждому VCF
    df = update_with_vcfs(df, args.vcfs)

    # Шаг 3: сохранить
    df.to_excel(args.output, index=False)  # движок openpyxl по умолчанию
    print(f"Готово: {args.output}")

if __name__ == "__main__":
    main()
