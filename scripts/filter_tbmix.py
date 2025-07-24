#!/usr/bin/env python3
"""
Фильтрация таблицы tbmix на основе шаблонов из tblg.

Использование:
    python filter_tbmix.py -f tbmix.tsv -t tblg.tsv -o result.tsv

Скрипт связывает строки по колонке `Sample`, оставляя в каждом из уровней
(Lvl1…Lvl5) только те теги, которые совпадают с шаблоном из tblg (символ
`*` соответствует любым символам). Дроби внутри скобок сохраняются.
"""
import argparse
import pandas as pd
import re
import sys
from pathlib import Path
from typing import Tuple, List, Optional

LEVELS = [f"Lvl{i}" for i in range(1, 6)]
ALT_LEVELS = {f"level_{i}": f"Lvl{i}" for i in range(1, 6)}  # tblg может иметь такие имена

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Фильтрует уровни Lvl1–Lvl5 из tbmix по шаблонам из tblg (поддержка '*')."
    )
    p.add_argument('-f', '--tbmix', required=True, help='TSV‑файл tbmix')
    p.add_argument('-t', '--tblg',  required=True, help='TSV‑файл tblg (шаблоны)')
    p.add_argument('-o', '--output', required=True, help='Куда сохранить результат (TSV)')
    return p.parse_args()

# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def read_table(path: str | Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep='\t', dtype=str).fillna('')
    except Exception as e:
        sys.exit(f"Ошибка при чтении '{path}': {e}")
    df.columns = df.columns.str.strip()
    return df


def write_table(df: pd.DataFrame, path: str | Path) -> None:
    try:
        df.to_csv(path, sep='\t', index=False)
    except Exception as e:
        sys.exit(f"Ошибка при записи '{path}': {e}")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def compile_regex(pattern: str) -> re.Pattern:
    """Создаёт регулярное выражение из шаблона с `*`."""
    return re.compile('^' + re.escape(pattern).replace(r'\*', '.*') + '$')


def parse_part(part: str) -> Optional[Tuple[str, str]]:
    """Разбивает `ТЕГ(дроби)` на (тег, дроби), поддерживая скобки внутри тега."""
    m = re.match(r'^(.*)\(([^()]*)\)$', part.strip())
    if m:
        tag = m.group(1).strip()
        frac = m.group(2)
        return tag, frac
    return None


def filter_level(mix_value: str, pattern: str) -> str:
    """Возвращает строку mix_value, содержащую только записи, совпавшие с pattern."""
    # Пустой или NaN шаблон → пустая строка
    if not pattern or pd.isna(pattern):
        return ''

    pattern = str(pattern)
    regex = compile_regex(pattern)
    kept: List[str] = []
    for part in (p for p in str(mix_value).split(';') if p.strip()):
        parsed = parse_part(part)
        if not parsed:
            continue
        tag, frac = parsed
        if regex.match(tag):
            kept.append(f"{pattern}({frac})")
    return ';'.join(kept)

# ---------------------------------------------------------------------------
# Core
# ---------------------------------------------------------------------------

def filter_mix(tbmix: pd.DataFrame, tblg: pd.DataFrame) -> pd.DataFrame:
    # Подготовка tblg: переименуем level_1 → Lvl1 и добавим недостающие колонки
    tblg = tblg.rename(columns=ALT_LEVELS)
    for lvl in LEVELS:
        if lvl not in tblg.columns:
            tblg[lvl] = ''

    # Left join по Sample
    merged = tbmix.merge(
        tblg[['Sample', *LEVELS]], on='Sample', how='left', suffixes=('', '_pat')
    ).fillna('')  # NaN → '' для простоты

    result = tbmix.copy()

    for lvl in LEVELS:
        pat_col = f"{lvl}_pat"
        result[lvl] = [
            filter_level(mix_val, pat_val)
            for mix_val, pat_val in zip(merged[lvl], merged[pat_col])
        ]

    return result

# ---------------------------------------------------------------------------
# Entry‑point
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    tbmix = read_table(args.tbmix)
    tblg  = read_table(args.tblg)

    filtered = filter_mix(tbmix, tblg)

    write_table(filtered, args.output)
    print(f"Готово! Фильтрованная таблица сохранена в '{args.output}'")


if __name__ == '__main__':
    main()