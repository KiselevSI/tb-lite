#!/usr/bin/env python3
"""
tb_summary.py – write TB‑Profiler results into a nice two‑level Excel sheet
without pandas' MultiIndex export bugs.

Usage examples
--------------
# На экран:
python tb_summary.py -i results/*/*.json

# В Excel:
python tb_summary.py -i results/*/*.json -o dr.xlsx
"""

from pathlib import Path
import json
import argparse
from typing import Dict, Any, List, Tuple

import pandas as pd
import xlsxwriter  # engine we’ll use directly

# ─────────── постоянные данные ───────────────────────────────────────────── #
DRUGS = [
    "Rifampicin", "Isoniazid", "Ethambutol", "Pyrazinamide",
    "Moxifloxacin", "Levofloxacin", "Bedaquiline", "Delamanid",
    "Pretomanid", "Linezolid", "Streptomycin", "Amikacin",
    "Kanamycin", "Capreomycin", "Clofazimine", "Ethionamide",
    "Para-aminosalicylic_acid", "Cycloserine",
]
FIELDS = ["Pos", "GeneName", "Mutation", "Freq", "Confidence"]
DRUG_MAP = {d.lower(): d for d in DRUGS}

# ─────────── helpers ─────────────────────────────────────────────────────── #
def load_json(path: Path) -> Dict[str, Any]:
    with path.open() as fh:
        return json.load(fh)

def summarise(sample: Dict[str, Any]) -> Dict[Tuple[str, str], str]:
    """One JSON → dict {(drug, field): value} with ';' accumulation."""
    row = {(d, f): "" for d in DRUGS for f in FIELDS}

    for var in sample.get("dr_variants", []):
        for drug in var.get("drugs", []):
            d = DRUG_MAP.get(drug["drug"].strip().lower())
            if not d:
                continue

            def acc(field: str, val):
                if not val:  # skip blanks
                    return
                k = (d, field)
                row[k] = f"{row[k]};{val}" if row[k] else str(val)

            acc("Pos", var.get("pos"))
            acc("GeneName", var.get("gene_name"))
            acc("Mutation", var.get("change"))
            acc("Freq", var.get("freq"))
            acc("Confidence", drug.get("confidence"))

    return row

def write_excel(df: pd.DataFrame, out_path: Path) -> None:
    """Write Excel so that 'Sample' is on the 2‑nd header level."""
    wb = xlsxwriter.Workbook(out_path)
    ws = wb.add_worksheet("TB‑Summary")

    hdr_fmt  = wb.add_format({"bold": True, "align": "center",
                              "valign": "vcenter", "border": 1})
    sub_fmt  = wb.add_format({"align": "center", "valign": "vcenter",
                              "border": 1})
    cell_fmt = wb.add_format({"border": 1})

    # ── header rows ──
    ws.write_blank(0, 0, "", hdr_fmt)      # верхний‑левый угол (пусто)
    ws.write(1, 0, "Sample", sub_fmt)      # «Sample» только на 2‑м уровне

    col = 1
    for drug in DRUGS:
        start, end = col, col + len(FIELDS) - 1
        ws.merge_range(0, start, 0, end, drug, hdr_fmt)  # 1‑й уровень
        for i, field in enumerate(FIELDS):               # 2‑й уровень
            ws.write(1, start + i, field, sub_fmt)
        col = end + 1

    # ── data rows ──
    for r, (sample_id, row) in enumerate(df.iterrows(), start=2):
        ws.write(r, 0, sample_id, cell_fmt)
        col = 1
        for drug in DRUGS:
            for field in FIELDS:
                ws.write(r, col, row[(drug, field)], cell_fmt)
                col += 1

    ws.freeze_panes(2, 1)
    ws.autofilter(1, 0, df.shape[0] + 1, len(DRUGS) * len(FIELDS))
    wb.close()

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="TB‑Profiler 2‑level Excel summary.")
    p.add_argument("-i", "--input", required=True, nargs="+", metavar="FILE.json",
                   help="One or more TB‑Profiler JSON files.")
    p.add_argument("-o", "--output", metavar="FILE.xlsx",
                   help="Output Excel file; leave blank to just print.")
    return p.parse_args()

# ─────────── main ────────────────────────────────────────────────────────── #
def main() -> None:
    args = parse_args()
    rows, ids = [], []

    for fname in args.input:
        p = Path(fname)
        sid = p.stem.split(".")[0]
        rows.append(summarise(load_json(p)))
        ids.append(sid)

    cols = pd.MultiIndex.from_product([DRUGS, FIELDS])
    df = pd.DataFrame(rows, index=ids, columns=cols)
    df.index.name = "Sample"  # keeps table tidy when printing

    if args.output:
        write_excel(df, Path(args.output))
        print(f"Saved {len(df)} sample(s) ➜ {args.output}")
    else:
        flat = df.copy()
        flat.columns = [".".join(c) for c in flat.columns]
        print(flat.reset_index())

if __name__ == "__main__":
    main()
