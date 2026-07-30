#!/usr/bin/env python3
"""Сборка нормализованных IS6110-таблиц из вывода ISMapper.

На вход подаются таблицы вызовов ISMapper ``<sample>__<ref>_table.txt`` и,
опционально, соседние файлы удалённых хитов ``<sample>__<ref>_table`` (без
расширения). Если файл удалённых хитов не передан явно, он ищется рядом с
таблицей вызовов на диске.

Формат ``_table.txt`` (TSV с шапкой)::

    region  orientation  x  y  gap  call  percent_ID  percent_cov
    left_gene  left_description  left_strand  left_distance
    right_gene right_description right_strand right_distance  gene_interruption

Формат ``_table`` (удалённые хиты, TSV с шапкой)::

    left_flank  right_flank  removal_reason  per_id  coverage  comparison_type

где ``left_flank`` -- диапазон вида ``3122027-3122203``.

На выходе -- шесть TSV, которые напрямую принимает импортёр TB Platform
``backend/scripts/import_is6110_tsv.py`` (порядок колонок и имена совпадают с
его словарём ``COLS``), плюс четыре диагностических файла:

* ``is_element.tsv``           -- справочник IS-элементов (одна строка: IS6110)
* ``is6110_site.tsv``          -- уникальные сайты вставки по всей когорте
* ``ismapper_run.tsv``         -- запуск ISMapper на образец/референс
* ``sample_is6110_site.tsv``   -- вызовы: образец -> сайт
* ``is6110_removed_hit.tsv``   -- отброшенные ISMapper кандидаты
* ``is6110_sample_summary.tsv``-- сводка по образцу
* ``import_stats.tsv``, ``missing_samples.tsv``, ``duplicate_runs.tsv``,
  ``is6110_import_warnings.tsv`` -- диагностика

Пустые значения записываются пустой строкой: COPY-импортёр трактует их как NULL.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

IS_ELEMENT_ID = 1
IS_ELEMENT_NAME = "IS6110"
TOOL_NAME = "ISMapper"
PARAMS_JSON = json.dumps({"element_name": IS_ELEMENT_NAME}, separators=(",", ":"))

# <sample>__<reference>_table.txt -- таблица вызовов.
# Файл без .txt рядом с ней -- удалённые хиты.
TABLE_RE = re.compile(r"^(?P<sample>.+?)__(?P<ref>.+)_table\.txt$")
REMOVED_RE = re.compile(r"^(?P<sample>.+?)__(?P<ref>.+)_table$")

COLUMNS: Dict[str, List[str]] = {
    "is_element": ["id", "name", "sequence_md5", "length_bp"],
    "is6110_site": [
        "id", "is_element_id", "reference_acc", "x_pos", "y_pos", "gap", "orientation",
        "left_gene", "left_description", "left_strand", "left_distance",
        "right_gene", "right_description", "right_strand", "right_distance",
        "gene_interruption",
    ],
    "ismapper_run": [
        "id", "sample_name", "is_element_id", "reference_acc", "tool_name", "tool_version",
        "params_json", "output_dir", "table_txt_path", "removed_hits_path",
        "total_calls", "confident_calls", "imprecise_calls", "uncertain_calls",
        "gene_interruptions",
    ],
    "sample_is6110_site": [
        "id", "run_id", "sample_name", "site_id", "region_name", "call_raw", "call_type",
        "is_imprecise", "is_uncertain", "percent_id", "percent_cov",
    ],
    "is6110_removed_hit": [
        "id", "run_id", "left_flank_start", "left_flank_end",
        "right_flank_start", "right_flank_end",
        "removal_reason", "percent_id", "coverage", "comparison_type",
    ],
    "is6110_sample_summary": [
        "sample_name", "run_id", "reference_acc", "total_calls", "confident_calls",
        "imprecise_calls", "uncertain_calls", "gene_interruptions", "table_txt_path",
    ],
}

# Поля аннотации сайта: сравниваются при повторной встрече того же сайта.
SITE_ANNOTATION_FIELDS = [
    "left_gene", "left_description", "left_strand", "left_distance",
    "right_gene", "right_description", "right_strand", "right_distance",
    "gene_interruption",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build normalized IS6110 tables from ISMapper output.",
    )
    parser.add_argument(
        "--input-list",
        help="Text file with one ISMapper table path per line.",
    )
    parser.add_argument(
        "-i",
        "--input",
        nargs="*",
        default=[],
        help="ISMapper table paths (_table.txt and optionally _table).",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Directory for the generated TSV files.",
    )
    return parser.parse_args()


def read_input_list(path: Path) -> List[Path]:
    items: List[Path] = []
    seen: set[str] = set()

    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            key = str(Path(line))
            if key in seen:
                continue

            seen.add(key)
            items.append(Path(line))

    return items


def merge_inputs(args: argparse.Namespace) -> List[Path]:
    items: List[Path] = []
    seen: set[str] = set()

    sources: List[Path] = [Path(item) for item in args.input]
    if args.input_list:
        sources.extend(read_input_list(Path(args.input_list)))

    for path in sources:
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        items.append(path)

    return items


def empty_to_blank(value: Optional[str]) -> str:
    """Пустое/отсутствующее значение -> пустая строка (NULL для COPY)."""
    if value is None:
        return ""
    return value.strip()


def parse_int(value: Optional[str]) -> Optional[int]:
    text = empty_to_blank(value)
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def parse_float(value: Optional[str]) -> Optional[float]:
    text = empty_to_blank(value)
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def format_number(value) -> str:
    return "" if value is None else str(value)


def format_bool(value: Optional[bool]) -> str:
    if value is None:
        return ""
    return "true" if value else "false"


def parse_bool(value: Optional[str]) -> Optional[bool]:
    text = empty_to_blank(value).lower()
    if not text:
        return None
    if text in ("true", "yes", "1"):
        return True
    if text in ("false", "no", "0"):
        return False
    return None


def call_flags(raw_value: Optional[str]) -> Tuple[str, str, bool, bool]:
    """'novel*' -> ('novel*', 'novel', imprecise=True, uncertain=False).

    Как в TB Platform (`backend/run_results_service.py::_is6110_call_flags`):
    '*' -- неточный вызов, '?' -- неопределённый.
    """
    raw = empty_to_blank(raw_value)
    if not raw:
        return "", "", False, False
    is_imprecise = raw.endswith("*")
    is_uncertain = raw.endswith("?")
    call_type = raw.rstrip("*?").strip()
    return raw, call_type, is_imprecise, is_uncertain


def parse_flank(value: Optional[str]) -> Tuple[Optional[int], Optional[int]]:
    """'3122027-3122203' -> (3122027, 3122203)."""
    text = empty_to_blank(value)
    if not text:
        return None, None
    parts = text.split("-")
    if len(parts) != 2:
        return None, None
    return parse_int(parts[0]), parse_int(parts[1])


WARNING_COLUMNS = ["level", "type", "sample_name", "reference_acc", "path", "message"]


class TsvWriter:
    """Потоковая запись TSV: когорта на 190k образцов даёт миллионы строк,
    держать их в памяти нельзя (процесс IS6110_TABLES ограничен 4 ГБ)."""

    def __init__(self, path: Path, header: List[str]) -> None:
        self.handle = path.open("w", encoding="utf-8", newline="")
        self.writer = csv.writer(self.handle, delimiter="\t", lineterminator="\n")
        self.writer.writerow(header)
        self.count = 0

    def write(self, row: List[str]) -> None:
        self.writer.writerow(row)
        self.count += 1

    def warn(self, kind: str, sample: str, reference: str, path: str, message: str) -> None:
        self.write(["WARNING", kind, sample, reference, path, message])

    def close(self) -> None:
        self.handle.close()


def discover_runs(
    inputs: List[Path],
) -> Tuple[List[Tuple[str, str, Path, Optional[Path]]], List[List[str]]]:
    """Собрать (sample, reference, table.txt, removed_hits) и список дублей.

    Дубль -- вторая и последующие таблицы для одной пары (образец, референс).
    Оставляем первую в лексикографическом порядке пути: результат детерминирован
    независимо от порядка стейджинга Nextflow.
    """
    # Файл удалённых хитов привязываем к каталогу: у одного образца может быть
    # несколько каталогов ISMapper (например, при повторном запуске).
    explicit_removed: Dict[Path, Path] = {}
    tables: Dict[Tuple[str, str], List[Path]] = {}

    for path in inputs:
        name = path.name
        match = TABLE_RE.match(name)
        if match:
            key = (match.group("sample"), match.group("ref"))
            tables.setdefault(key, []).append(path)
            continue

        if REMOVED_RE.match(name):
            explicit_removed.setdefault(path.resolve(), path)

    runs: List[Tuple[str, str, Path, Optional[Path]]] = []
    duplicates: List[List[str]] = []

    for key in sorted(tables):
        sample, reference = key
        candidates = sorted(tables[key], key=lambda item: str(item))
        kept = candidates[0]
        for skipped in candidates[1:]:
            duplicates.append([
                sample,
                reference,
                str(kept.resolve()),
                str(skipped.resolve()),
                "duplicate (sample_name, reference_acc)",
            ])

        sibling = kept.parent / kept.name[: -len(".txt")]
        removed = explicit_removed.get(sibling.resolve())
        if removed is None and sibling.exists():
            removed = sibling

        runs.append((sample, reference, kept, removed))

    return runs, duplicates


def read_table_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_tables(
    runs: List[Tuple[str, str, Path, Optional[Path]]],
    writers: Dict[str, TsvWriter],
    warnings: TsvWriter,
) -> None:
    """Разобрать все прогоны, записывая строки в TSV по мере чтения.

    В памяти остаётся только индекс сайтов: ключ -> id и хеш аннотации для
    обнаружения расхождений. Сами строки никуда не накапливаются.
    """
    site_ids: Dict[Tuple, int] = {}
    site_annotation_hashes: Dict[int, int] = {}

    call_id = 0
    removed_id = 0

    for run_id, (sample, reference, table_path, removed_path) in enumerate(runs, start=1):
        rows = read_table_rows(table_path)

        total = 0
        confident = 0
        imprecise = 0
        uncertain = 0
        interruptions = 0

        # line_number: 1 -- шапка, поэтому данные начинаются со второй строки.
        for offset, row in enumerate(rows, start=2):
            orientation = empty_to_blank(row.get("orientation"))
            x_pos = parse_int(row.get("x"))
            y_pos = parse_int(row.get("y"))
            if orientation not in ("F", "R") or x_pos is None or y_pos is None:
                warnings.warn(
                    "bad_main_row",
                    sample,
                    reference,
                    str(table_path.resolve()),
                    f"Skipped line {offset}: missing/invalid orientation/x/y",
                )
                continue

            gap = parse_int(row.get("gap"))
            annotation = {
                "left_gene": empty_to_blank(row.get("left_gene")),
                "left_description": empty_to_blank(row.get("left_description")),
                "left_strand": format_number(parse_int(row.get("left_strand"))),
                "left_distance": format_number(parse_int(row.get("left_distance"))),
                "right_gene": empty_to_blank(row.get("right_gene")),
                "right_description": empty_to_blank(row.get("right_description")),
                "right_strand": format_number(parse_int(row.get("right_strand"))),
                "right_distance": format_number(parse_int(row.get("right_distance"))),
                "gene_interruption": format_bool(parse_bool(row.get("gene_interruption"))),
            }

            annotation_values = [annotation[field] for field in SITE_ANNOTATION_FIELDS]
            annotation_hash = hash(tuple(annotation_values))

            site_key = (IS_ELEMENT_ID, reference, x_pos, y_pos, gap, orientation)
            site_id = site_ids.get(site_key)
            if site_id is None:
                site_id = len(site_ids) + 1
                site_ids[site_key] = site_id
                site_annotation_hashes[site_id] = annotation_hash
                writers["is6110_site"].write([
                    str(site_id),
                    str(IS_ELEMENT_ID),
                    reference,
                    str(x_pos),
                    str(y_pos),
                    format_number(gap),
                    orientation,
                ] + annotation_values)
            elif site_annotation_hashes[site_id] != annotation_hash:
                warnings.warn(
                    "site_annotation_conflict",
                    sample,
                    reference,
                    str(table_path.resolve()),
                    f"Site id {site_id} has conflicting annotation; kept first occurrence",
                )

            raw_call, call_type, row_imprecise, row_uncertain = call_flags(row.get("call"))

            call_id += 1
            writers["sample_is6110_site"].write([
                str(call_id),
                str(run_id),
                sample,
                str(site_id),
                empty_to_blank(row.get("region")),
                raw_call,
                call_type,
                format_bool(row_imprecise),
                format_bool(row_uncertain),
                format_number(parse_float(row.get("percent_ID"))),
                format_number(parse_float(row.get("percent_cov"))),
            ])

            total += 1
            if row_uncertain:
                uncertain += 1
            elif row_imprecise:
                imprecise += 1
            else:
                confident += 1
            if annotation["gene_interruption"] == "true":
                interruptions += 1

        if removed_path is not None:
            for row in read_table_rows(removed_path):
                left_start, left_end = parse_flank(row.get("left_flank"))
                right_start, right_end = parse_flank(row.get("right_flank"))
                removed_id += 1
                writers["is6110_removed_hit"].write([
                    str(removed_id),
                    str(run_id),
                    format_number(left_start),
                    format_number(left_end),
                    format_number(right_start),
                    format_number(right_end),
                    empty_to_blank(row.get("removal_reason")),
                    format_number(parse_float(row.get("per_id"))),
                    format_number(parse_float(row.get("coverage"))),
                    empty_to_blank(row.get("comparison_type")),
                ])

        output_dir = str(table_path.resolve().parent)
        table_txt_path = str(table_path.resolve())
        removed_hits_path = str(removed_path.resolve()) if removed_path is not None else ""

        writers["ismapper_run"].write([
            str(run_id),
            sample,
            str(IS_ELEMENT_ID),
            reference,
            TOOL_NAME,
            "",
            PARAMS_JSON,
            output_dir,
            table_txt_path,
            removed_hits_path,
            str(total),
            str(confident),
            str(imprecise),
            str(uncertain),
            str(interruptions),
        ])

        writers["is6110_sample_summary"].write([
            sample,
            str(run_id),
            reference,
            str(total),
            str(confident),
            str(imprecise),
            str(uncertain),
            str(interruptions),
            table_txt_path,
        ])


def write_tsv(path: Path, header: List[str], rows: List[List[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    inputs = merge_inputs(args)
    runs, duplicates = discover_runs(inputs)

    writers = {
        table: TsvWriter(outdir / f"{table}.tsv", header)
        for table, header in COLUMNS.items()
    }
    warnings = TsvWriter(outdir / "is6110_import_warnings.tsv", WARNING_COLUMNS)
    try:
        writers["is_element"].write([str(IS_ELEMENT_ID), IS_ELEMENT_NAME, "", ""])
        build_tables(runs, writers, warnings)
    finally:
        for writer in writers.values():
            writer.close()
        warnings.close()

    write_tsv(
        outdir / "duplicate_runs.tsv",
        ["sample_name", "reference_acc", "kept_table_txt_path", "skipped_table_txt_path", "reason"],
        duplicates,
    )
    write_tsv(outdir / "missing_samples.tsv", ["sample_name"], [])

    samples = {run[0] for run in runs}
    stats = [
        ["found_samples_in_input_dirs", str(len(samples))],
        ["runs_selected", str(len(runs))],
        ["duplicate_runs_skipped", str(len(duplicates))],
        ["unique_is6110_sites", str(writers["is6110_site"].count)],
        ["sample_is6110_calls", str(writers["sample_is6110_site"].count)],
        ["removed_hits", str(writers["is6110_removed_hit"].count)],
        ["warnings", str(warnings.count)],
    ]
    write_tsv(outdir / "import_stats.tsv", ["metric", "value"], stats)


if __name__ == "__main__":
    main()
