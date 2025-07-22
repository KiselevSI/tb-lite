#!/usr/bin/env python3
"""
Create a *sample‑centric* wide table: **one row per readcounts file**, with a contiguous
block of columns for every variant defined in a single *levels.tsv* template.  Each
variant block repeats the fields

    Lineage, Level, Pos, Ref, Alt,
    Coverage, avg_bq, avg_mapq,
    Alt_fraction, Alt/Ref,
    A_fraction, C_fraction, G_fraction, T_fraction

Column names are suffixed with the 1‑based variant index (e.g. `Pos_1`, `Alt_fraction_3`).
The very first column is always `Sample` (stem of the readcounts filename).

Typical usage:

```bash
python tb-lineage.py \
    -t reference_levels.tsv \
    -i sample1_readcounts.tsv sample2_readcounts.tsv sample3_readcounts.tsv \
    [-f 0.5]
```

Arguments:
    -t, --table   one *levels* TSV (columns POS/REF/ALT or Pos/Ref/Alt)
    -i, --inputs  one or more *readcounts* TSV files (pos/ref/depth + base counts)
    -f, --filter  optional float; drop rows where **no** variant's Alt_fraction
                  exceeds threshold (keeps row if *any* variant > threshold)

The script prints the resulting TSV to *stdout*.
"""

import argparse
from pathlib import Path
import sys
import pandas as pd

# -----------------------------------------------------------------------------
# CLI parsing
# -----------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Merge one levels.tsv with multiple readcounts TSVs into a per‑sample wide table."
    )
    p.add_argument("-t", "--table", required=True, metavar="LEVELS.tsv")
    p.add_argument("-i", "--inputs", required=True, nargs="+", metavar="READCOUNTS.tsv")
    p.add_argument("-f", "--filter", type=float, default=None, metavar="THRESHOLD")
    return p.parse_args()

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

BASES = ("A", "C", "G", "T", "N")
VAR_FIELDS = [
    "Lineage",
    "Level",
    "Pos",
    "Ref",
    "Alt",
    "Coverage",
    "avg_bq",
    "avg_mapq",
    "Alt_fraction",
    "Alt/Ref",
    "A_fraction",
    "C_fraction",
    "G_fraction",
    "T_fraction",
]


def load_readcounts(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = df.rename(columns={"pos": "Pos", "ref": "Ref", "depth": "Coverage"})
    for col in ["Coverage", "avg_bq", "avg_mapq", *BASES]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def load_levels(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = df.rename(
        columns={
            "POS": "Pos",
            "REF": "Ref",
            "ALT": "Alt",
            "lineage": "Lineage",
            "LEVEL": "Level",
            "level": "Level",
        }
    )
    return df


def compute_metrics(level_row: pd.Series, rc_subset: pd.Series) -> dict:
    """Given one variant row from levels and matching readcounts row, compute metrics."""
    alt = str(level_row["Alt"]).upper()
    ref = str(level_row["Ref"]).upper()
    cov = rc_subset["Coverage"]

    base_counts = {b: rc_subset.get(b) for b in BASES}
    alt_cnt = base_counts.get(alt)
    ref_cnt = base_counts.get(ref)

    def frac(cnt):
        return pd.NA if pd.isna(cov) or cov == 0 or pd.isna(cnt) else round(cnt / cov, 2)

    metrics = {
        "Lineage": level_row["Lineage"],
        "Level": level_row["Level"],
        "Pos": level_row["Pos"],
        "Ref": level_row["Ref"],
        "Alt": level_row["Alt"],
        "Coverage": cov,
        "avg_bq": rc_subset["avg_bq"],
        "avg_mapq": rc_subset["avg_mapq"],
        "Alt_fraction": frac(alt_cnt),
        "Alt/Ref": pd.NA if pd.isna(ref_cnt) or ref_cnt == 0 or pd.isna(alt_cnt) else round(alt_cnt / ref_cnt, 2),
        "A_fraction": frac(base_counts.get("A")),
        "C_fraction": frac(base_counts.get("C")),
        "G_fraction": frac(base_counts.get("G")),
        "T_fraction": frac(base_counts.get("T")),
    }
    return metrics

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    args = parse_args()
    levels_path = Path(args.table)
    if not levels_path.exists():
        sys.stderr.write(f"Error: levels file '{levels_path}' not found.\n")
        sys.exit(1)

    levels_df = load_levels(levels_path)
    # Assign a stable variant index (1‑based) to keep column ordering consistent
    levels_df = levels_df.reset_index(drop=True)
    levels_df["VarIdx"] = levels_df.index + 1

    sample_rows = []
    for rc_path in [Path(p) for p in args.inputs]:
        if not rc_path.exists():
            sys.stderr.write(f"Warning: readcounts file '{rc_path}' not found, skipping.\n")
            continue
        rc_df = load_readcounts(rc_path)
        sample_name = rc_path.stem

        # Merge readcounts with levels to align rows; left join keeps variant order
        merged = pd.merge(levels_df, rc_df, on=["Pos", "Ref"], how="left", validate="1:1")

        row_data = {"Sample": sample_name}
        for _, row in merged.iterrows():
            idx = row["VarIdx"]
            metrics = compute_metrics(row, row)
            for key, val in metrics.items():
                row_data[f"{key}_{idx}"] = val
        sample_rows.append(row_data)

    if not sample_rows:
        sys.stderr.write("Error: no valid readcounts files processed.\n")
        sys.exit(1)

    wide_df = pd.DataFrame(sample_rows)

    # Optional filter: keep sample if any Alt_fraction_i > threshold
    if args.filter is not None:
        alt_cols = [c for c in wide_df.columns if c.startswith("Alt_fraction_")]
        if alt_cols:
            wide_df = wide_df[wide_df[alt_cols].gt(args.filter).any(axis=1)]

    # Sort columns: Sample first, then blocks by variant idx in ascending order
    block_cols = []
    for idx in range(1, len(levels_df) + 1):
        block_cols.extend([f"{field}_{idx}" for field in VAR_FIELDS])
    output_cols = ["Sample", *block_cols]
    wide_df = wide_df.reindex(columns=output_cols)
    for col in wide_df.columns:
        # подходит: колонка содержит "fraction" или начинается с "Alt/Ref"
        if "fraction" in col or col.startswith("Alt/Ref"):
            wide_df[col] = pd.to_numeric(wide_df[col], errors="coerce")

    # --- Write to stdout with comma as decimal separator ---
    wide_df.to_csv(sys.stdout, sep="\t", index=False, decimal=",")


if __name__ == "__main__":
    main()