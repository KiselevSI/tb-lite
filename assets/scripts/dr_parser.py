#!/usr/bin/env python3
"""
Updated tb_summary.py – write TB‑Profiler results into a nice two‑level Excel sheet
without pandas' MultiIndex export bugs.

Key changes in this version
---------------------------
1. **Freq rounding** – all frequency values are rounded to **two decimal places** so that
   numbers such as `67.33333333` become `67.33`.
2. **Other Variants column** – a new column that shows the **count of objects in the
   `other_variants` list** inside each TB‑Profiler JSON.  The column appears just after
   the *Sample* column in both the console output and the Excel file, and it stays
   frozen together with Sample when you scroll.

Usage examples remain the same:

```bash
# To stdout:
python tb_summary.py -i results/*/*.json

# To Excel:
python tb_summary.py -i results/*/*.json -o dr.xlsx
```
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
OTHER_VARIANTS_COL = ("", "Other Variants")  # MultiIndex tuple

# ─────────── helpers ─────────────────────────────────────────────────────── #

def load_json(path: Path) -> Dict[str, Any]:
    with path.open() as fh:
        return json.load(fh)


def summarise(sample: Dict[str, Any]) -> Dict[Tuple[str, str], str]:
    """Convert one TB‑Profiler JSON → dict suitable for a MultiIndex DataFrame.

    • Drug‑specific data are stored under the key ``dr_variants``.
    • "Other Variants" is simply the length of the top‑level
      ``other_variants`` list (0 if absent).
    • Frequency values are rounded to **two decimal places**.
    """
    # Start with empty strings for every drug/field pair
    row: Dict[Tuple[str, str], str] = {
        (d, f): "" for d in DRUGS for f in FIELDS
    }

    # Pre‑populate the special column
    row[OTHER_VARIANTS_COL] = str(len(sample.get("other_variants", [])))

    # Helper to accumulate semicolon‑separated strings
    def acc(key: Tuple[str, str], val: Any):
        if val in ("", None):  # skip blanks
            return
        current = row[key]
        row[key] = f"{current};{val}" if current else str(val)

    for var in sample.get("dr_variants", []):
        # Iterate over all drugs declared for this variant
        for drug in var.get("drugs", []):
            d = DRUG_MAP.get(drug["drug"].strip().lower())
            if not d:
                continue  # skip unrecognised drugs

            # Extract frequency with rounding
            freq_val = var.get("freq")
            if freq_val not in (None, ""):
                try:
                    freq_val = f"{float(freq_val):.2f}"
                except (ValueError, TypeError):
                    # Leave as‑is if conversion fails
                    pass

            # Accumulate all fields
            acc((d, "Pos"), var.get("pos"))
            acc((d, "GeneName"), var.get("gene_name"))
            acc((d, "Mutation"), var.get("change"))
            acc((d, "Freq"), freq_val)
            acc((d, "Confidence"), drug.get("confidence"))

    return row


def write_excel(df: pd.DataFrame, out_path: Path) -> None:
    """Write the summary to Excel with a two‑level header.

    The first column is *Sample*, the second is *Other Variants*, followed by
    a block for every drug and its five per‑variant fields.
    """
    wb = xlsxwriter.Workbook(out_path)
    ws = wb.add_worksheet("TB‑Summary")

    hdr_fmt = wb.add_format({"bold": True, "align": "center", "valign": "vcenter", "border": 1})
    sub_fmt = wb.add_format({"align": "center", "valign": "vcenter", "border": 1})
    cell_fmt = wb.add_format({"border": 1})

    # ── header rows ──
    ws.write_blank(0, 0, "", hdr_fmt)        # (0,0) — top‑left corner
    ws.write(1, 0, "Sample", sub_fmt)         # row 1 header for Sample

    ws.write_blank(0, 1, "", hdr_fmt)        # blank top‑level for Other Variants
    ws.write(1, 1, "Other Variants", sub_fmt)

    col = 2  # start writing drug blocks from the third column
    for drug in DRUGS:
        start, end = col, col + len(FIELDS) - 1
        ws.merge_range(0, start, 0, end, drug, hdr_fmt)  # top‑level header
        for i, field in enumerate(FIELDS):  # second‑level headers
            ws.write(1, start + i, field, sub_fmt)
        col = end + 1

    # ── data rows ──
    for r, (sample_id, row) in enumerate(df.iterrows(), start=2):
        ws.write(r, 0, sample_id, cell_fmt)
        ws.write(r, 1, row[OTHER_VARIANTS_COL], cell_fmt)
        col = 2
        for drug in DRUGS:
            for field in FIELDS:
                ws.write(r, col, row[(drug, field)], cell_fmt)
                col += 1

    # Keep Sample and Other Variants columns visible when scrolling
    ws.freeze_panes(2, 2)
    # Apply autofilter across the full width of the data
    ws.autofilter(1, 0, df.shape[0] + 1, df.shape[1])
    wb.close()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="TB‑Profiler 2‑level Excel summary.")
    p.add_argument("-i", "--input", required=True, nargs="+", metavar="FILE.json", help="One or more TB‑Profiler JSON files.")
    p.add_argument("-o", "--output", metavar="FILE.xlsx", help="Output Excel file; leave blank to just print.")
    return p.parse_args()

# ─────────── main ────────────────────────────────────────────────────────── #

def main() -> None:
    args = parse_args()
    rows: List[Dict[Tuple[str, str], str]] = []
    ids: List[str] = []

    for fname in args.input:
        p = Path(fname)
        sid = p.stem.split(".")[0]
        rows.append(summarise(load_json(p)))
        ids.append(sid)

    # Build MultiIndex columns: first the extra column, then the drug/field grid
    columns = pd.MultiIndex.from_tuples([OTHER_VARIANTS_COL]).append(
        pd.MultiIndex.from_product([DRUGS, FIELDS])
    )

    df = pd.DataFrame(rows, index=ids, columns=columns)
    df.index.name = "Sample"

    if args.output:
        write_excel(df, Path(args.output))
        print(f"Saved {len(df)} sample(s) ➜ {args.output}")
    else:
        flat = df.copy()
        # Flatten MultiIndex for pretty console output
        flat.columns = [c[1] if not c[0] else "{}.{}".format(*c) for c in flat.columns]
        print(flat.reset_index())


if __name__ == "__main__":
    main()
