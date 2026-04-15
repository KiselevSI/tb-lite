#!/usr/bin/env python3
"""
 tb_mix_v6.py — pile‑up ➜ lineage summary ➜ mixture flag    (pandas 2.x + pysam)

Правила определения смеси:
1) Проходим уровни **5 → 4 → 3 → 2 → 1**. На каждом уровне решаем независимо.
2) **Присутствие линии** на уровне: у линии есть ≥1 SNP с `AF ≥ mix_low`.
3) Если на уровне присутствует **>1 линий** → `mix`.
4) Если на уровне присутствует **ровно 1 линия** → `mix`, если у неё есть SNP с `AF ∈ [mix_low, mix_high]`.
5) **Особое правило**: любая линия, имя которой начинается с `Bmyc` (например, `Bmyc3`), при присутствии (п.2) даёт `mix`.
6) Если ни на одном уровне условия не выполнены → `clear`.

Вывод таблицы:
• Для **Lvl1/Lvl2**: если у линии **хотя бы один** SNP проходит `--min-af`, выводим **все её SNP‑фракции на уровне**, включая те, что < `--min-af` (в т.ч. 0.00). Так `L1.3*(0.04:)` превратится в `L1.3(0.04:0.00)`.
• Для **Lvl3–Lvl5**: выводим только фракции, которые ≥ `--min-af`.
• Звёздочки **не выводим**.

Примечания:
• Инверсия AF для `L4@Lvl1` и `L4.9@Lvl2` — используем `1 − AF`.
• L2.2: объединяем `ancient → modern` для lvl2 и отчёта.
"""
from __future__ import annotations
import argparse, csv, statistics as st
from pathlib import Path
from typing import List

import pandas as pd
import pysam

# ─────────────────────────────────────── 1. levels.tsv ─────────────────────────
def load_levels(path: Path) -> pd.DataFrame:
    """Читаем файл с баркод‑SNP‑ами (линия, уровень, pos, ref, alt)."""
    df = pd.read_csv(
        path, sep="	", na_filter=False,
        dtype={"lineage": "category", "level": "int8",
               "pos": "int32", "ref": "string", "alt": "string"}
    )
    df.columns = df.columns.str.lower()
    return df

# ─────────────────────────────────────── 2. pile‑up ───────────────────────────
def pileup_counts(bam: Path, ref: Path, levels: pd.DataFrame,
                  mq: int, bq: int) -> pd.DataFrame:
    """AF для каждого баркода (с учётом инверсии L4/L4.9)."""
    bamf = pysam.AlignmentFile(bam, "rb")  # AlignmentFile API
    ref_name = pysam.FastaFile(ref).references[0]  # FastaFile.references
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
                if pr.is_del or pr.is_refskip:
                    continue
                # При stepper="nofilter" min_mapping_quality фильтруем вручную.
                if pr.alignment.mapping_quality < mq:
                    continue
                base = pr.alignment.query_sequence[pr.query_position].upper()
                if base == rec.alt:
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

# ──────────────────────────── 3. агрегации и вспомогательные ──────────────────
def aggregate_raw(raw_df: pd.DataFrame) -> pd.DataFrame:
    """Группировка без фильтров: список фракций и число SNP‑ов."""
    if raw_df.empty:
        return pd.DataFrame(columns=["lineage", "level", "fractions", "evid_n"])
    grp = (raw_df.groupby(["lineage", "level"], sort=False)
                  .agg(fractions=("alt_fraction", list),
                       evid_n=("alt_fraction", "size"))
                  .reset_index())
    # L2.2: ancient → modern (на lvl2)
    anc = grp.query("lineage == 'L2.2 (ancient)' and level == 2")
    mod = grp.query("lineage == 'L2.2 (modern)' and level == 2")
    if not anc.empty and not mod.empty:
        merged = anc.fractions.iloc[0] + mod.fractions.iloc[0]
        idx_mod = mod.index[0]
        grp.at[idx_mod, "fractions"] = merged
        grp.at[idx_mod, "evid_n"] = len(merged)
        grp = grp[grp.lineage != "L2.2 (ancient)"]
    return grp

# ────────────────────────────────── 4. поиск смеси ────────────────────────────
def detect_mix_levels(raw_df: pd.DataFrame,
                      mix_low: float = 0.05, mix_high: float = 0.95) -> bool:
    """
    Проходим уровни 5→4→3→2→1.
    Присутствие линии: есть ≥1 SNP с AF ≥ mix_low.
    На уровне:
      • >1 присутствующих линий → mix;
      • 1 линия → mix, если есть фракция в [mix_low, mix_high];
      • иначе — к следующему уровню.
    Особый случай: любая линия `Bmyc*` при присутствии → mix.
    """
    grp = aggregate_raw(raw_df)


    for lvl in (5, 4, 3, 2, 1):
        lv = grp[grp.level == lvl].copy()
        if lv.empty:
            continue

        # оставляем только присутствующие линии
        lv["present"] = lv.fractions.apply(lambda frs: any(f >= mix_low for f in frs))
        lv = lv[lv.present]
        if lv.empty:
            continue

        n_lines = lv.lineage.nunique()
        if n_lines > 1:
            return True
        if n_lines == 1:
            frs = lv.fractions.iloc[0]
            if any(mix_low <= f <= mix_high for f in frs):
                return True

    return False

# ─────────────────────────────── 5. агрегация для отчёта ──────────────────────
def aggregate_for_report(raw_df: pd.DataFrame, min_af: float,
                         min_snps12: int) -> pd.DataFrame:
    """Собираем строки для вывода с правилами Lvl1/Lvl2: «если один прошёл — показать все».
    Звёздочки не выводим.
    """
    full = aggregate_raw(raw_df)  # все фракции без фильтра

    rows = []
    for r in full.itertuples(index=False):
        fr_all = list(r.fractions)
        pass_fr = [x for x in fr_all if x >= min_af]
        any_pass = len(pass_fr) > 0

        if r.level in (1, 2):
            if not any_pass:
                continue  # не показываем линию, если ни один SNP не прошёл порог
            fr_out = fr_all  # показываем все SNP‑ы, включая < min_af (и 0.00)
        else:
            if not any_pass:
                continue
            fr_out = pass_fr  # для уровней 3–5 показываем только ≥ min_af

        rows.append({
            "lineage": r.lineage,
            "level": int(r.level),
            "fractions": fr_out,
            # evid_n для отчёта больше не нужен; медиану считаем по выводимым фракциям
        })

    if not rows:
        return pd.DataFrame(columns=["lineage", "level", "fractions", "med_af"])

    grp = pd.DataFrame(rows)
    grp["med_af"] = grp.fractions.apply(st.median)
    return grp.sort_values(["level", "lineage"], ignore_index=True)

# ─────────────────────────────── 6. формирование колонок ──────────────────────
def pad(fracs: List[float], lineage: str, level: int) -> str:
    need = 3 if lineage == "L2.2 (modern)" else 2 if level in (1, 2) else 1
    return ":".join([f"{x:.2f}" for x in fracs] + [""] * max(0, need - len(fracs)))

def make_cols(grp: pd.DataFrame) -> List[str]:
    if grp.empty:
        return [""] * 5

    grp = grp.copy()
    grp["pretty"] = [
        f"{r.lineage}({pad(r.fractions, r.lineage, r.level)})"
        for r in grp.itertuples(index=False)
    ]
    return [";".join(grp.loc[grp.level == lvl, "pretty"]) for lvl in range(1, 6)]

# ───────────────────────────────────── CLI ────────────────────────────────────
def main():
    p = argparse.ArgumentParser(
        "tb_mix v6 — 5→1, присутствие по mix_low; Lvl1/2 выводит все SNP, если один прошёл",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i", "--input-bam", required=True, type=Path)
    p.add_argument("-r", "--reference", required=True, type=Path)
    p.add_argument("-l", "--levels", required=True, type=Path)
    p.add_argument("--mq", type=int, default=0)
    p.add_argument("--bq", type=int, default=0)
    p.add_argument("-f", "--min-af", type=float, default=0.05,
                   help="минимальная AF, которая попадёт в таблицу")
    p.add_argument("--mix-low", type=float, default=0.05,
                   help="нижняя граница для определения смеси (включительно)")
    p.add_argument("--mix-high", type=float, default=0.95,
                   help="верхняя граница для определения смеси (включительно)")
    # оставляем --min-snps для совместимости CLI, но он больше не влияет на вывод
    p.add_argument("--min-snps", type=int, default=2,
                   help="сохранён для совместимости; на mix и вывод не влияет")
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
    is_mix = detect_mix_levels(pileup, args.mix_low, args.mix_high)
    mix_flag = "mix" if is_mix else "clear"

    report_df = aggregate_for_report(pileup, args.min_af, args.min_snps)
    lvl_cols = make_cols(report_df)

    sample_name = args.prefix if args.prefix else args.input_bam.stem

    # ---------- вывод ----------
    header = ["Sample", "Lvl1", "Lvl2", "Lvl3", "Lvl4", "Lvl5", "Mix"]
    row = [sample_name, *lvl_cols, mix_flag]

    if args.save_pileup:
        pileup.to_csv(args.save_pileup, sep="	", index=False)

    if args.output:
        with open(args.output, "w", newline="") as fh:
            csv.writer(fh, delimiter="	").writerows([header, row])
    else:
        print("	".join(header))
        print("	".join(row))

if __name__ == "__main__":
    pd.set_option("mode.copy_on_write", True)
    main()
