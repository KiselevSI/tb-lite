#!/usr/bin/env python3
"""
tb_mix_v3.py — pile‑up ➜ lineage summary ➜ mixture flag    (pandas 2.x + pysam)

• Микс = ≥2 линий на одном уровне, median(AF) ∈ [0.10, 0.90].
• ‑f/--min‑af влияет только на отображение, а не на определение смеси.
"""
from __future__ import annotations
import argparse, csv, statistics as st
from pathlib import Path
from typing import List  # PEP 484 typing 🔗 :contentReference[oaicite:1]{index=1}

import pandas as pd        # pandas groupby/median 🔗 :contentReference[oaicite:2]{index=2}
import pysam               # pile‑up API 🔗 :contentReference[oaicite:3]{index=3}

# ─────────────────────────────────────── 1. levels.tsv ─────────────────────────
def load_levels(path: Path) -> pd.DataFrame:
    """Читаем файл с баркод‑SNP‑ами (линия, уровень, pos, ref, alt)."""
    df = pd.read_csv(
        path, sep="\t", na_filter=False,
        dtype={"lineage": "category", "level": "int8",
               "pos": "int32", "ref": "string", "alt": "string"}
    )
    df.columns = df.columns.str.lower()
    return df

# ─────────────────────────────────────── 2. pile‑up ───────────────────────────
def pileup_counts(bam: Path, ref: Path, levels: pd.DataFrame,
                  mq: int, bq: int) -> pd.DataFrame:
    """AF для каждого баркода (с учётом инверсии L4/L4.9)."""
    bamf = pysam.AlignmentFile(bam, "rb")
    ref_name = pysam.FastaFile(ref).references[0]
    rows: List[dict] = []

    for rec in levels.itertuples(index=False):
        depth = alt = 0
        for col in bamf.pileup(
            ref_name, rec.pos - 1, rec.pos,
            truncate=True, min_base_quality=bq, stepper="nofilter"
        ):
            if col.reference_pos != rec.pos - 1:
                continue
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip or pr.alignment.mapping_quality < mq:
                    continue
                if pr.alignment.query_sequence[pr.query_position].upper() == rec.alt:
                    alt += 1
                depth += 1
        af = alt / depth if depth else 0.0

        # инверсия долей для L4@Lvl1 и L4.9@Lvl2
        if (rec.lineage == "L4" and rec.level == 1) or \
           (rec.lineage == "L4.9" and rec.level == 2):
            af = 1 - af

        rows.append({"lineage": rec.lineage, "level": rec.level,
                     "alt_fraction": af})
    bamf.close()
    return pd.DataFrame(rows)

# ────────────────────────────────── 3. поиск смеси ────────────────────────────
def detect_mix(raw_df: pd.DataFrame, min_snps12: int,
               low: float = 0.10, high: float = 0.90) -> bool:
    """
    Микс, если где‑то ≥2 линии и median(AF) в диапазоне [low, high].
    Смотрим на все SNP‑ы, без фильтра ‑f/--min-af.
    """
    grp = (raw_df.groupby(["lineage", "level"], sort=False)
                   .agg(med_af=("alt_fraction", "median"),
                        evid_n=("alt_fraction", "size"))
                   .reset_index())

    # учитываем только линии с достаточным числом баркодов и правильным AF
    mask = (grp.med_af.between(low, high)) & (grp.evid_n >= min_snps12)
    eligible = grp[mask]

    # на каком‑то уровне появилось ≥2 различных линии?
    by_level = eligible.groupby("level").lineage.nunique()  # nunique 🔗 :contentReference[oaicite:4]{index=4}
    return (by_level >= 2).any()

# ─────────────────────────────── 4. агрегация для отчёта ──────────────────────
def aggregate_for_report(df: pd.DataFrame, min_af: float,
                         min_snps12: int) -> pd.DataFrame:
    """Отфильтровать слабые SNP‑ы для табличного вывода."""
    df = df[df.alt_fraction >= min_af]
    if df.empty:
        return pd.DataFrame(columns=["lineage", "level", "fractions",
                                     "evid_n", "star", "med_af"])

    grp = (df.groupby(["lineage", "level"], sort=False)
             .agg(fractions=("alt_fraction", list),
                  evid_n=("alt_fraction", "size"))
             .reset_index())

    grp["star"] = grp.apply(
        lambda r: "*" if r.level in (1, 2) and 0 < r.evid_n < min_snps12 else "",
        axis=1)
    grp["med_af"] = grp.fractions.apply(st.median)          # statistics.median 🔗 :contentReference[oaicite:5]{index=5}

    # «ancient → modern» для L2.2
    anc = grp.query("lineage == 'L2.2 (ancient)'")
    mod = grp.query("lineage == 'L2.2 (modern)'")
    if (not anc.empty and not mod.empty
            and anc.evid_n.iloc[0] >= 2 and mod.evid_n.iloc[0] == 1):
        merged = anc.fractions.iloc[0] + mod.fractions.iloc[0]
        idx_mod = mod.index[0]
        grp.at[idx_mod, "fractions"] = merged
        grp.at[idx_mod, "star"] = ""
        grp.at[idx_mod, "med_af"] = st.median(merged)
        grp = grp[grp.lineage != "L2.2 (ancient)"]

    return grp.sort_values(["level", "lineage"], ignore_index=True)

# ─────────────────────────────── 5. формирование колонок ──────────────────────
def pad(fracs: List[float], lineage: str, level: int) -> str:
    need = 3 if lineage == "L2.2 (modern)" else 2 if level in (1, 2) else 1
    return ":".join([f"{x:.2f}" for x in fracs] + [""] * max(0, need - len(fracs)))

def make_cols(grp: pd.DataFrame) -> List[str]:
    grp = grp.assign(
        pretty=lambda g: g.apply(
            lambda r: f"{r.lineage}{r.star}({pad(r.fractions, r.lineage, r.level)})", axis=1))
    # пять уровней (пустые строки, если линии нет)
    return [";".join(grp.loc[grp.level == lvl, "pretty"]) for lvl in range(1, 6)]

# ───────────────────────────────────── CLI ────────────────────────────────────
def main():
    p = argparse.ArgumentParser("tb_mix v3 — отдельные пороги для mix и вывода")
    p.add_argument("-i", "--input-bam", required=True, type=Path)
    p.add_argument("-r", "--reference", required=True, type=Path)
    p.add_argument("-l", "--levels", required=True, type=Path)
    p.add_argument("--mq", type=int, default=0)
    p.add_argument("--bq", type=int, default=0)
    p.add_argument("-f", "--min-af", type=float, default=0.05,
                   help="минимальная AF, которая попадёт в таблицу (по умолчанию 0.05)")
    p.add_argument("--min-snps", type=int, default=2,
                   help="минимум SNP‑ов уровня 1/2 для учёта линии")
    p.add_argument("-s", "--save-pileup", type=Path,
                   help="сохранять полный pile‑up (без фильтра) в TSV")
    p.add_argument("-o", "--output", type=Path,
                   help="куда писать отчёт (таб‑разделённый); stdout, если не указано")
    p.add_argument("-p", "--prefix", type=str,
                   help="явное имя образца для первой колонки")
    args = p.parse_args()

    levels = load_levels(args.levels)
    pileup = pileup_counts(args.input_bam, args.reference, levels, args.mq, args.bq)

    # ---------- mix / report ----------
    mix_flag = "mix" if detect_mix(pileup, args.min_snps) else "clear"
    report_df = aggregate_for_report(pileup, args.min_af, args.min_snps)
    lvl_cols = make_cols(report_df)

    sample_name = args.prefix if args.prefix else args.input_bam.stem

    # ---------- вывод ----------
    header = ["Sample", "Lvl1", "Lvl2", "Lvl3", "Lvl4", "Lvl5", "Mix"]
    row = [sample_name, *lvl_cols, mix_flag]

    if args.save_pileup:
        pileup.to_csv(args.save_pileup, sep="\t", index=False)

    if args.output:
        with open(args.output, "w", newline="") as fh:
            csv.writer(fh, delimiter="\t").writerows([header, row])
    else:
        print("\t".join(header))
        print("\t".join(row))

if __name__ == "__main__":
    pd.set_option("mode.copy_on_write", True)   # Copy‑on‑write 🔗 :contentReference[oaicite:7]{index=7}
    main()
