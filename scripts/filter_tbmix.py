#!/usr/bin/env python3
"""
Фильтрация таблицы tbmix на основе шаблонов из tblg
с учётом эквивалентности суффиксов “(ancient)” и “(modern)”.

Использование:
    python filter_tbmix.py -f tbmix.tsv -t tblg.tsv -o result.tsv

Скрипт связывает строки по колонке `Sample`, оставляя в каждом из уровней
(Lvl1…Lvl5) только те теги, которые совпадают с шаблоном из tblg.
Поддержка:

* Символ `*` соответствует любым символам (как в shell‑маске).
* Если шаблон содержит тег `Lx.y (ancient)` и в tbmix встречается
  `Lx.y (modern)` (или наоборот), запись считается совпавшей и
  в выходной файл попадает в «ancient»‑варианте (то есть строго по шаблону).
* Дроби внутри скобок (пример: `0.77:0.81`) всегда сохраняются.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys
from typing import List, Optional, Tuple

import pandas as pd


# ---------------------------------------------------------------------------
# Константы
# ---------------------------------------------------------------------------

LEVELS = [f"Lvl{i}" for i in range(1, 6)]
ALT_LEVELS = {f"level_{i}": f"Lvl{i}" for i in range(1, 6)}  # tblg может иметь такие имена

# регулярка для распознавания суффиксов эпох
_EPOCH_RE = re.compile(r"\s*\((?:ancient|modern)\)\s*$", re.I)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Фильтрует уровни Lvl1–Lvl5 из tbmix по шаблонам из tblg "
            "(поддержка '*', а также эквивалентность '(ancient)'/'(modern)')."
        )
    )
    parser.add_argument("-f", "--tbmix", required=True, help="TSV‑файл tbmix")
    parser.add_argument("-t", "--tblg", required=True, help="TSV‑файл tblg (шаблоны)")
    parser.add_argument("-o", "--output", required=True, help="Куда сохранить результат (TSV)")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def read_table(path: str | Path) -> pd.DataFrame:
    """Чтение TSV с заменой NaN на пустые строки."""
    try:
        df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    except Exception as exc:
        sys.exit(f"Ошибка при чтении '{path}': {exc}")
    df.columns = df.columns.str.strip()
    return df


def write_table(df: pd.DataFrame, path: str | Path) -> None:
    """Сохранение TSV без индекса."""
    try:
        df.to_csv(path, sep="\t", index=False)
    except Exception as exc:
        sys.exit(f"Ошибка при записи '{path}': {exc}")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def compile_regex(pattern: str) -> re.Pattern[str]:
    """Создаёт регулярное выражение из шаблона с `*`."""
    return re.compile("^" + re.escape(pattern).replace(r"\*", ".*") + "$")


def parse_part(part: str) -> Optional[Tuple[str, str]]:
    """
    Разбивает `ТЕГ(дроби)` → (тег, дроби).

    Допускает скобки внутри тега.
    """
    match = re.match(r"^(.*)\(([^()]*)\)$", part.strip())
    if match:
        tag = match.group(1).strip()
        frac = match.group(2)
        return tag, frac
    return None


def _strip_epoch(tag: str) -> str:
    """Убирает суффикс (ancient)/(modern) из тега."""
    return _EPOCH_RE.sub("", tag).strip()


def filter_level(mix_value: str, pattern: str) -> str:
    """
    Возвращает строку `mix_value`, содержащую только записи, совпавшие с `pattern`.

    Правила совпадения:
    * Полная проверка regexp на основе шаблона (с `*`).
    * Дополнительно: если базовые части тега (без суффикса эпохи)
      совпадают, запись также считается совпавшей.
    В выходную строку всегда попадает тег из `pattern` – чтобы, к примеру,
    `(modern)` преобразовалось в `(ancient)`, если такое сказано в tblg.
    """
    if not pattern or pd.isna(pattern):
        return ""

    pattern = str(pattern).strip()
    regex = compile_regex(pattern)
    base_pat = _strip_epoch(pattern)

    kept: List[str] = []
    for part in (p for p in str(mix_value).split(";") if p.strip()):
        parsed = parse_part(part)
        if not parsed:
            continue
        tag, frac = parsed

        # Совпадение по полному шаблону или по базовой части без (ancient)/(modern)
        if regex.match(tag) or _strip_epoch(tag) == base_pat:
            kept.append(f"{pattern}({frac})")

    return ";".join(kept)


# ---------------------------------------------------------------------------
# Core
# ---------------------------------------------------------------------------


def filter_mix(tbmix: pd.DataFrame, tblg: pd.DataFrame) -> pd.DataFrame:
    """
    Применяет фильтрацию к каждому уровню (Lvl1…Lvl5) по шаблонам из tblg.
    """
    # 1. Приводим названия колонок tblg к LvlX и добавляем недостающие.
    tblg = tblg.rename(columns=ALT_LEVELS)
    for lvl in LEVELS:
        if lvl not in tblg.columns:
            tblg[lvl] = ""

    # 2. left‑join по Sample, чтобы у каждой строки tbmix был набор шаблонов.
    merged = tbmix.merge(
        tblg[["Sample", *LEVELS]],
        on="Sample",
        how="left",
        suffixes=("", "_pat"),
    ).fillna("")

    # 3. Копируем tbmix и заменяем значения колонок.
    result = tbmix.copy()

    for lvl in LEVELS:
        pat_col = f"{lvl}_pat"
        result[lvl] = [
            filter_level(mix_val, pat_val) for mix_val, pat_val in zip(merged[lvl], merged[pat_col])
        ]

    return result


# ---------------------------------------------------------------------------
# Entry‑point
# ---------------------------------------------------------------------------


def main() -> None:
    args = parse_args()

    tbmix = read_table(args.tbmix)
    tblg = read_table(args.tblg)

    filtered = filter_mix(tbmix, tblg)

    write_table(filtered, args.output)
    print(f"Готово! Фильтрованная таблица сохранена в '{args.output}'")


if __name__ == "__main__":
    main()
