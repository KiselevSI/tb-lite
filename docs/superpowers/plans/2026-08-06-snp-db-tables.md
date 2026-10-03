# SNP-таблицы для TB Platform — план реализации

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Пайплайн tb-lite должен после прогона отдавать готовые к ручному импорту файлы для таблиц `snp_sites`, `sample_snp_alleles`, `vcf_table` и SQL-скрипт, который грузит их в БД TB Platform и достраивает `sample_snp_profiles`.

**Architecture:** Конвертер VCF→шарды запускается один раз на батч внутри `main.nf` (процесс `SNP_DB_SHARDS`), финальная склейка шардов — один раз в `batch_reports.nf` (процесс `SNP_DB_TABLES`). Пайплайн к БД не подключается: `site_id` раздаёт Postgres при импорте через staging по `site_key`.

**Tech Stack:** Nextflow DSL2, Python 3.12 (только stdlib), bash/coreutils/gawk, PostgreSQL 17, тесты — `unittest` из stdlib.

Спека: `docs/superpowers/specs/2026-08-06-snp-db-tables-design.md`.

## Global Constraints

- Формула `site_key` менять **нельзя**: sha1-hex от `"\x1f".join([chrom, pos, ref, alt, effect, impact, gene, gene_id, product_name, strand, feature_id, biotype, rank, hgvs_c, hgvs_p, cds_pos, aa_pos])`. Проверочные значения: pos 1977 → `7229c7508aa7fb285706a91c398ca4b22bdccee9`, pos 4013 → `f03d041b126f9db3c71d4f62fb698395201f8ee5`.
- Никакой фильтрации при извлечении: ни по типу варианта, ни по FILTER, ни по DP/AF, ни по маскам. Берём каждую ALT-аллель каждой записи VCF.
- `product_name` и `strand` всегда пустые (в per-sample VCF нет INFO `name`/`strand`). Заполнять из feature table запрещено — сломает `site_key`.
- Пустое значение в TSV = SQL NULL.
- Порядок колонок фиксирован:
  - sites: `site_key, chrom, pos, ref, alt, qual, effect, impact, gene, gene_id, product_name, strand, feature_id, biotype, rank, hgvs_c, hgvs_p, cds_pos, aa_pos`
  - alleles: `sample_id, site_key, allele, dp, af, qual, filter`
  - vcf_rows: `sample_id, pos, alt`
- Шарды пишутся **без строки заголовка**; заголовок появляется только в итоговых файлах.
- Только stdlib Python — конвертер работает в `tb-lite/tb-platform-tables:1.1` (`python:3.12-slim`), пересборка контейнера не нужна.
- Тесты запускаются из корня репозитория: `python3 -m unittest discover -s tests -t . -v`.
- Исходник для переноса логики: `/home/zerg/git/tb-platform/new_temp_data/vcf_to_snp_db_tsv_parallel_scan.py` (далее «исходный скрипт»).

## File Structure

| Файл | Ответственность |
|---|---|
| `bin/vcf_to_snp_shards.py` | извлечение из VCF + CLI, пишет шарды одного батча |
| `bin/merge_snp_shards.sh` | склейка шардов в итоговые файлы + манифест |
| `bin/write_import_snp_sql.py` | генерация `import_snp.sql` |
| `modules/local/snp_db/shards/main.nf` | процесс `SNP_DB_SHARDS` |
| `modules/local/snp_db/tables/main.nf` | процесс `SNP_DB_TABLES` |
| `subworkflows/local/snp_db.nf` | обвязка `SNP_DB` |
| `tests/data/ERR4797591.mini.ann.vcf` | фикстура: 7 записей из реального прогона |
| `tests/test_vcf_to_snp_shards.py` | юнит-тесты конвертера |
| `tests/test_merge_snp_shards.py` | тесты склейки |
| `tests/test_import_snp_sql.py` | интеграционный тест SQL против локальной БД |
| `tests/nf/snp_db_shards_smoke.nf` | smoke-workflow для процесса `SNP_DB_SHARDS` |
| `tests/nf/snp_db_final_smoke.nf` | smoke-workflow для `SNP_DB_FINAL` (ветка пересборки шардов) |

---

### Task 1: Ядро конвертера — извлечение из VCF

**Files:**
- Create: `tests/data/ERR4797591.mini.ann.vcf`
- Create: `bin/vcf_to_snp_shards.py`
- Test: `tests/test_vcf_to_snp_shards.py`

**Interfaces:**
- Consumes: ничего.
- Produces: `extract_from_vcf(path: pathlib.Path) -> ExtractResult`, где `ExtractResult` — `NamedTuple(sample_id: str, sites: list[list[str]], alleles: list[list[str]], vcf_rows: list[list[str]])`. Порядок колонок в каждом списке — как в Global Constraints. Также модуль экспортирует константы `SITE_COLUMNS`, `ALLELE_COLUMNS`, `VCF_COLUMNS` (списки имён колонок).

- [ ] **Step 1: Собрать фикстуру из реального прогона**

Записи выбраны так, чтобы покрыть все типы вариантов и запись с двумя ANN-аннотациями (перекрывающиеся гены pks6/Rv0406c):

```bash
cd /home/zerg/git/tb-lite
mkdir -p tests/data
RUN=/mnt/data1/tb-lite-runs/bc1286ed-5e2e-42db-a06c-0a5a02258824/results/annotate_vcf/ERR4797591/ERR4797591.annotated.ann.vcf
{
  grep '^#' "$RUN"
  awk -F'\t' '$2==1977 || $2==4013 || $2==26747 || $2==26957 || $2==125830 || $2==489935 || $2==3133054' "$RUN"
} > tests/data/ERR4797591.mini.ann.vcf
grep -vc '^#' tests/data/ERR4797591.mini.ann.vcf   # ожидается: 7
```

Состав: 1977 (snp, intergenic), 4013 (snp, missense), 26747 (del), 26957 (complex), 125830 (ins), 489935 (snp, 2 ANN), 3133054 (mnp).

- [ ] **Step 2: Написать падающий тест**

`tests/test_vcf_to_snp_shards.py`:

```python
import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "bin"))

from vcf_to_snp_shards import ALLELE_COLUMNS, SITE_COLUMNS, VCF_COLUMNS, extract_from_vcf

FIXTURE = REPO_ROOT / "tests" / "data" / "ERR4797591.mini.ann.vcf"

# Значения из продовой таблицы snp_sites (tb_database.dump). Совместимость по
# site_key обязательна: по нему новые данные схлопываются с существующими.
GOLDEN_SITE_KEYS = {
    1977: "7229c7508aa7fb285706a91c398ca4b22bdccee9",
    4013: "f03d041b126f9db3c71d4f62fb698395201f8ee5",
}


class ExtractFromVcfTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = extract_from_vcf(FIXTURE)

    def test_sample_id_from_vcf_header(self):
        self.assertEqual(self.result.sample_id, "ERR4797591")

    def test_site_key_matches_production(self):
        by_pos = {int(row[2]): row[0] for row in self.result.sites}
        for pos, expected in GOLDEN_SITE_KEYS.items():
            self.assertEqual(by_pos[pos], expected, f"site_key разошёлся на pos {pos}")

    def test_column_counts(self):
        self.assertTrue(all(len(row) == len(SITE_COLUMNS) for row in self.result.sites))
        self.assertTrue(all(len(row) == len(ALLELE_COLUMNS) for row in self.result.alleles))
        self.assertTrue(all(len(row) == len(VCF_COLUMNS) for row in self.result.vcf_rows))

    def test_keeps_every_variant_type(self):
        # 7 записей VCF, каждая с одной ALT -> 7 строк vcf_table
        self.assertEqual(len(self.result.vcf_rows), 7)
        self.assertEqual({int(row[1]) for row in self.result.vcf_rows},
                         {1977, 4013, 26747, 26957, 125830, 489935, 3133054})

    def test_multi_ann_record_expands_sites_but_not_vcf_rows(self):
        # pos 489935 перекрывается двумя генами -> 2 site/allele, но 1 строка vcf_table
        sites_489935 = [row for row in self.result.sites if row[2] == "489935"]
        alleles_489935 = [row for row in self.result.alleles
                          if row[1] in {row2[0] for row2 in sites_489935}]
        vcf_489935 = [row for row in self.result.vcf_rows if row[1] == "489935"]
        self.assertEqual(len(sites_489935), 2)
        self.assertEqual(len(alleles_489935), 2)
        self.assertEqual(len(vcf_489935), 1)
        self.assertEqual(len(self.result.sites), 8)
        self.assertEqual(len(self.result.alleles), 8)

    def test_product_name_and_strand_always_empty(self):
        pn = SITE_COLUMNS.index("product_name")
        strand = SITE_COLUMNS.index("strand")
        self.assertTrue(all(row[pn] == "" for row in self.result.sites))
        self.assertTrue(all(row[strand] == "" for row in self.result.sites))

    def test_allele_metrics_match_production_semantics(self):
        af = ALLELE_COLUMNS.index("af")
        filt = ALLELE_COLUMNS.index("filter")
        dp = ALLELE_COLUMNS.index("dp")
        self.assertTrue(all(row[af] == "1" for row in self.result.alleles))
        self.assertTrue(all(row[filt] == "" for row in self.result.alleles))
        self.assertTrue(all(row[dp].isdigit() for row in self.result.alleles))

    def test_site_keys_unique_within_sample(self):
        keys = [row[0] for row in self.result.sites]
        self.assertEqual(len(keys), len(set(keys)))


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 3: Запустить тест и убедиться, что он падает**

Run: `cd /home/zerg/git/tb-lite && python3 -m unittest tests.test_vcf_to_snp_shards -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'vcf_to_snp_shards'`

- [ ] **Step 4: Перенести ядро парсинга из исходного скрипта**

Создать `bin/vcf_to_snp_shards.py`. Скопировать **дословно** строки 109–413 исходного скрипта — это `ANN_SPLIT_RE`, `MISSING`, `eprint`, `clean`, `as_float`, `as_int`, `float_to_tsv`, `int_to_tsv`, `open_text`, `open_out`, `strip_vcf_suffix`, `parse_info`, `split_list_value`, `select_list_value`, `parse_format`, `gt_alt_indices`, `estimate_af_from_depths`, `extract_metrics`, `is_snp`, `ann_entries_for_alt`, `ann_field`, `site_row_from_ann`.

Из скопированного удалить `is_snp` (строки 292–294) — фильтрация по типу варианта запрещена. Остальное не трогать: любая правка `site_row_from_ann` ломает совместимость по `site_key`.

Шапка файла и импорты:

```python
#!/usr/bin/env python3
"""Конвертация per-sample аннотированных VCF в шарды для таблиц TB Platform.

Пишет три шарда без заголовка:
  sites.<tag>.tsv.gz     site_key + атрибуты снипа (site_id раздаёт Postgres)
  alleles.<tag>.tsv.gz   образец x site_key
  vcfrows.<tag>.tsv.gz   образец x позиция x альтернативная аллель

Фильтрации нет: берём каждую ALT-аллель каждой записи VCF.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import gzip
import hashlib
import os
import re
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Sequence
```

Константы колонок:

```python
SITE_COLUMNS = [
    "site_key", "chrom", "pos", "ref", "alt", "qual", "effect", "impact", "gene",
    "gene_id", "product_name", "strand", "feature_id", "biotype", "rank",
    "hgvs_c", "hgvs_p", "cds_pos", "aa_pos",
]
ALLELE_COLUMNS = ["sample_id", "site_key", "allele", "dp", "af", "qual", "filter"]
VCF_COLUMNS = ["sample_id", "pos", "alt"]

_QUAL_IDX = SITE_COLUMNS.index("qual")
```

- [ ] **Step 5: Реализовать `extract_from_vcf`**

Дописать в `bin/vcf_to_snp_shards.py`:

```python
class ExtractResult(NamedTuple):
    sample_id: str
    sites: List[List[str]]
    alleles: List[List[str]]
    vcf_rows: List[List[str]]


def _site_qual(row: Sequence[str]) -> float:
    value = as_float(row[_QUAL_IDX])
    return value if value is not None else float("-inf")


def extract_from_vcf(path: Path) -> ExtractResult:
    """Разобрать один per-sample аннотированный VCF.

    site_key считается по той же формуле, что и в проде, поэтому строки
    схлопываются с существующей таблицей snp_sites через ON CONFLICT.
    """
    sample_id = strip_vcf_suffix(path)
    sample_col_idx: Optional[int] = None

    sites: Dict[str, List[str]] = {}
    alleles: List[List[str]] = []
    seen_pairs: set = set()
    vcf_rows: List[List[str]] = []
    seen_vcf: set = set()

    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n\r").split("\t")
                sample_columns = header[9:] if len(header) > 9 else []
                if len(sample_columns) > 1:
                    raise ValueError(
                        f"Multi-sample VCF is not supported: {path} has "
                        f"{len(sample_columns)} sample columns"
                    )
                if sample_columns:
                    sample_id = sample_columns[0]
                    sample_col_idx = 9
                continue
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n\r").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _vid, ref, alt_text, qual_text, filt, info_text = parts[:8]

            alt_list = alt_text.split(",") if alt_text and alt_text != "." else []
            if not alt_list:
                continue

            qual = as_float(qual_text)
            info = parse_info(info_text)
            fmt_map: Dict[str, str] = {}
            if sample_col_idx is not None and len(parts) > sample_col_idx:
                fmt_map = parse_format(parts[8], parts[sample_col_idx])

            gt_indices = gt_alt_indices(fmt_map.get("GT", ""), len(alt_list))
            alt_indices = list(range(len(alt_list))) if gt_indices is None else gt_indices

            for alt_idx in alt_indices:
                alt = clean(alt_list[alt_idx])
                if not alt:
                    continue

                dp, af = extract_metrics(fmt_map, info, alt_idx)

                vcf_key = (clean(pos), alt)
                if vcf_key not in seen_vcf:
                    seen_vcf.add(vcf_key)
                    vcf_rows.append([sample_id, clean(pos), alt])

                single_alt = len(alt_list) == 1
                for ann_parts in ann_entries_for_alt(info, alt, single_alt=single_alt):
                    site_row = site_row_from_ann(chrom, pos, ref, alt, qual, info, ann_parts)
                    site_key = site_row[0]
                    previous = sites.get(site_key)
                    if previous is None or _site_qual(site_row) > _site_qual(previous):
                        sites[site_key] = site_row

                    pair = (sample_id, site_key)
                    if pair in seen_pairs:
                        continue
                    seen_pairs.add(pair)
                    alleles.append([
                        sample_id,
                        site_key,
                        alt,
                        int_to_tsv(dp),
                        float_to_tsv(af),
                        float_to_tsv(qual),
                        clean(filt),
                    ])

    return ExtractResult(sample_id, list(sites.values()), alleles, vcf_rows)
```

- [ ] **Step 6: Запустить тесты**

Run: `cd /home/zerg/git/tb-lite && python3 -m unittest tests.test_vcf_to_snp_shards -v`
Expected: PASS, 8 тестов

- [ ] **Step 7: Коммит**

```bash
cd /home/zerg/git/tb-lite
git add bin/vcf_to_snp_shards.py tests/test_vcf_to_snp_shards.py tests/data/ERR4797591.mini.ann.vcf
git commit -m "feat(snp-db): извлечение snp_sites/alleles/vcf_table из аннотированного VCF"
```

---

### Task 2: CLI конвертера — запись шардов батча

**Files:**
- Modify: `bin/vcf_to_snp_shards.py`
- Test: `tests/test_vcf_to_snp_shards.py`

**Interfaces:**
- Consumes: `extract_from_vcf`, `ExtractResult`, `SITE_COLUMNS`, `ALLELE_COLUMNS`, `VCF_COLUMNS` из Task 1.
- Produces: CLI `vcf_to_snp_shards.py --file-list <txt> --tag <tag> -o <dir> [-j N]`. В `<dir>` создаются `sites.<tag>.tsv.gz`, `alleles.<tag>.tsv.gz`, `vcfrows.<tag>.tsv.gz` (все без заголовка), `samples.<tag>.txt`, `counts.<tag>.tsv`. Формат `counts.<tag>.tsv` — две колонки без заголовка: `metric<TAB>value`, метрики `sites`, `alleles`, `vcfrows`, `samples`.

- [ ] **Step 1: Написать падающий тест**

Дописать в `tests/test_vcf_to_snp_shards.py`:

```python
import gzip
import subprocess
import tempfile


def _read_gz_rows(path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return [row for row in csv.reader(handle, delimiter="\t") if row]


class ShardCliTest(unittest.TestCase):
    """Два образца с одинаковыми сайтами, но разным QUAL: проверяем дедупликацию."""

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        tmp = Path(cls._tmp.name)

        # Второй образец: то же содержимое, другое имя и вдвое меньший QUAL.
        second = tmp / "SAMPLE2.annotated.ann.vcf"
        text = FIXTURE.read_text(encoding="utf-8")
        lines = []
        for line in text.splitlines():
            if line.startswith("#CHROM"):
                lines.append(line.replace("ERR4797591", "SAMPLE2"))
            elif line.startswith("#"):
                lines.append(line)
            else:
                parts = line.split("\t")
                parts[5] = f"{float(parts[5]) / 2:.2f}"
                lines.append("\t".join(parts))
        second.write_text("\n".join(lines) + "\n", encoding="utf-8")

        file_list = tmp / "vcfs.txt"
        file_list.write_text(f"{FIXTURE}\n{second}\n", encoding="utf-8")

        cls.outdir = tmp / "out"
        subprocess.run(
            [sys.executable, str(REPO_ROOT / "bin" / "vcf_to_snp_shards.py"),
             "--file-list", str(file_list), "--tag", "batch_1",
             "-o", str(cls.outdir), "-j", "2"],
            check=True,
        )

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_writes_expected_files(self):
        for name in ("sites.batch_1.tsv.gz", "alleles.batch_1.tsv.gz",
                     "vcfrows.batch_1.tsv.gz", "samples.batch_1.txt",
                     "counts.batch_1.tsv"):
            self.assertTrue((self.outdir / name).exists(), name)

    def test_shards_have_no_header(self):
        first = _read_gz_rows(self.outdir / "alleles.batch_1.tsv.gz")[0]
        self.assertNotEqual(first[0], "sample_id")

    def test_row_counts(self):
        # 8 уникальных site_key на два образца, по 8 аллелей и 7 vcf-строк на образец
        self.assertEqual(len(_read_gz_rows(self.outdir / "sites.batch_1.tsv.gz")), 8)
        self.assertEqual(len(_read_gz_rows(self.outdir / "alleles.batch_1.tsv.gz")), 16)
        self.assertEqual(len(_read_gz_rows(self.outdir / "vcfrows.batch_1.tsv.gz")), 14)

    def test_sites_keep_max_qual(self):
        rows = _read_gz_rows(self.outdir / "sites.batch_1.tsv.gz")
        by_pos = {int(row[2]): float(row[5]) for row in rows}
        # У ERR4797591 QUAL выше вдвое, значит должен победить он
        self.assertAlmostEqual(by_pos[1977], 1643.47, places=2)

    def test_samples_file(self):
        samples = (self.outdir / "samples.batch_1.txt").read_text().split()
        self.assertEqual(sorted(samples), ["ERR4797591", "SAMPLE2"])

    def test_counts_file(self):
        counts = dict(
            line.split("\t")
            for line in (self.outdir / "counts.batch_1.tsv").read_text().splitlines()
        )
        self.assertEqual(counts["sites"], "8")
        self.assertEqual(counts["alleles"], "16")
        self.assertEqual(counts["vcfrows"], "14")
        self.assertEqual(counts["samples"], "2")
```

Добавить `import csv` в начало файла тестов.

- [ ] **Step 2: Запустить тест и убедиться, что он падает**

Run: `cd /home/zerg/git/tb-lite && python3 -m unittest tests.test_vcf_to_snp_shards.ShardCliTest -v`
Expected: FAIL — `subprocess.CalledProcessError` (у скрипта ещё нет CLI)

- [ ] **Step 3: Реализовать воркер и сборку шардов**

Дописать в `bin/vcf_to_snp_shards.py`:

```python
def _write_rows(path: Path, rows: Sequence[Sequence[str]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


def _worker(task):
    index, path_str, tmpdir_str = task
    tmpdir = Path(tmpdir_str)
    try:
        result = extract_from_vcf(Path(path_str))
    except Exception as exc:  # шард одного образца не должен ронять весь батч
        return {"index": index, "sample_id": None, "error": f"{type(exc).__name__}: {exc}"}
    _write_rows(tmpdir / f"sites_{index:09d}.tsv", result.sites)
    _write_rows(tmpdir / f"alleles_{index:09d}.tsv", result.alleles)
    _write_rows(tmpdir / f"vcfrows_{index:09d}.tsv", result.vcf_rows)
    return {"index": index, "sample_id": result.sample_id, "error": None}


def _read_paths(file_list: Path) -> List[Path]:
    paths: List[Path] = []
    seen: set = set()
    with open(file_list, "r", encoding="utf-8") as handle:
        for line in handle:
            raw = line.strip()
            if not raw or raw.startswith("#"):
                continue
            key = os.path.abspath(raw)
            if key in seen:
                continue
            seen.add(key)
            paths.append(Path(key))
    return paths


def _merge_shards(tmpdir: Path, outdir: Path, tag: str) -> Dict[str, int]:
    sites: Dict[str, List[str]] = {}
    for shard in sorted(tmpdir.glob("sites_*.tsv")):
        with open(shard, "r", encoding="utf-8", newline="") as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if not row:
                    continue
                previous = sites.get(row[0])
                if previous is None or _site_qual(row) > _site_qual(previous):
                    sites[row[0]] = row

    counts = {"sites": 0, "alleles": 0, "vcfrows": 0}
    with gzip.open(outdir / f"sites.{tag}.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for row in sites.values():
            writer.writerow(row)
            counts["sites"] += 1

    for kind in ("alleles", "vcfrows"):
        with gzip.open(outdir / f"{kind}.{tag}.tsv.gz", "wt", encoding="utf-8", newline="") as out:
            for shard in sorted(tmpdir.glob(f"{kind}_*.tsv")):
                with open(shard, "r", encoding="utf-8") as handle:
                    for line in handle:
                        if line.strip():
                            out.write(line)
                            counts[kind] += 1
    return counts


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--file-list", required=True, type=Path,
                        help="Текстовый файл со списком путей к VCF, по одному на строку.")
    parser.add_argument("--tag", required=True, help="Тег шарда, попадает в имена файлов.")
    parser.add_argument("-o", "--out", required=True, type=Path, help="Каталог вывода.")
    parser.add_argument("-j", "--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    args = parser.parse_args(argv)

    paths = _read_paths(args.file_list)
    if not paths:
        eprint(f"ERROR: пустой список VCF: {args.file_list}")
        return 1

    args.out.mkdir(parents=True, exist_ok=True)
    samples: List[str] = []
    errors: List[str] = []

    with tempfile.TemporaryDirectory(dir=str(args.out)) as tmpdir_str:
        tmpdir = Path(tmpdir_str)
        tasks = [(i, str(path), tmpdir_str) for i, path in enumerate(paths)]
        with cf.ProcessPoolExecutor(max_workers=args.jobs) as pool:
            for outcome in pool.map(_worker, tasks, chunksize=8):
                if outcome["error"]:
                    errors.append(f"{paths[outcome['index']]}: {outcome['error']}")
                else:
                    samples.append(outcome["sample_id"])
        counts = _merge_shards(tmpdir, args.out, args.tag)

    samples = sorted(set(samples))
    (args.out / f"samples.{args.tag}.txt").write_text(
        "".join(f"{name}\n" for name in samples), encoding="utf-8"
    )
    counts["samples"] = len(samples)
    (args.out / f"counts.{args.tag}.tsv").write_text(
        "".join(f"{key}\t{value}\n" for key, value in counts.items()), encoding="utf-8"
    )

    for message in errors:
        eprint(f"WARNING: {message}")
    eprint(
        f"[{args.tag}] samples={counts['samples']} sites={counts['sites']} "
        f"alleles={counts['alleles']} vcfrows={counts['vcfrows']} errors={len(errors)}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Запустить все тесты конвертера**

Run: `cd /home/zerg/git/tb-lite && python3 -m unittest tests.test_vcf_to_snp_shards -v`
Expected: PASS, 14 тестов

- [ ] **Step 5: Коммит**

```bash
cd /home/zerg/git/tb-lite
git add bin/vcf_to_snp_shards.py tests/test_vcf_to_snp_shards.py
git commit -m "feat(snp-db): CLI конвертера, шарды батча с дедупликацией по site_key"
```

---

### Task 3: Склейка шардов в итоговые таблицы

**Files:**
- Create: `bin/merge_snp_shards.sh`
- Test: `tests/test_merge_snp_shards.py`

**Interfaces:**
- Consumes: шарды из Task 2 (`sites.*.tsv.gz`, `alleles.*.tsv.gz`, `vcfrows.*.tsv.gz`, `samples.*.txt`, `counts.*.tsv`).
- Produces: `merge_snp_shards.sh --shard-dir <dir> --out-dir <dir> [--threads N] [--sort-mem 2G]`. В `<out-dir>` появляются `snp_sites.tsv.gz`, `sample_snp_alleles.tsv.gz`, `vcf_table.tsv.gz` (все **с** заголовком), `samples.txt`, `snp_db_manifest.tsv`. Манифест — две колонки без заголовка `metric<TAB>value`, метрики: `shards`, `samples`, `sites_rows_before_dedup`, `sites_rows`, `alleles_rows`, `vcf_rows`.

- [ ] **Step 1: Написать падающий тест**

`tests/test_merge_snp_shards.py`:

```python
import csv
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "bin" / "merge_snp_shards.sh"
CONVERTER = REPO_ROOT / "bin" / "vcf_to_snp_shards.py"
FIXTURE = REPO_ROOT / "tests" / "data" / "ERR4797591.mini.ann.vcf"


def _read_gz_rows(path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return [row for row in csv.reader(handle, delimiter="\t") if row]


class MergeSnpShardsTest(unittest.TestCase):
    """Два шарда с пересекающимися сайтами и разным QUAL."""

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        tmp = Path(cls._tmp.name)
        shard_dir = tmp / "shards"

        second = tmp / "SAMPLE2.annotated.ann.vcf"
        lines = []
        for line in FIXTURE.read_text(encoding="utf-8").splitlines():
            if line.startswith("#CHROM"):
                lines.append(line.replace("ERR4797591", "SAMPLE2"))
            elif line.startswith("#"):
                lines.append(line)
            else:
                parts = line.split("\t")
                parts[5] = f"{float(parts[5]) * 3:.2f}"
                lines.append("\t".join(parts))
        second.write_text("\n".join(lines) + "\n", encoding="utf-8")

        for tag, vcf in (("batch_1", FIXTURE), ("batch_2", second)):
            listing = tmp / f"{tag}.txt"
            listing.write_text(f"{vcf}\n", encoding="utf-8")
            subprocess.run(
                [sys.executable, str(CONVERTER), "--file-list", str(listing),
                 "--tag", tag, "-o", str(shard_dir), "-j", "1"],
                check=True,
            )

        cls.outdir = tmp / "out"
        subprocess.run(
            ["bash", str(SCRIPT), "--shard-dir", str(shard_dir),
             "--out-dir", str(cls.outdir), "--threads", "2"],
            check=True,
        )

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_outputs_exist(self):
        for name in ("snp_sites.tsv.gz", "sample_snp_alleles.tsv.gz",
                     "vcf_table.tsv.gz", "samples.txt", "snp_db_manifest.tsv"):
            self.assertTrue((self.outdir / name).exists(), name)

    def test_headers_present(self):
        self.assertEqual(_read_gz_rows(self.outdir / "snp_sites.tsv.gz")[0][0], "site_key")
        self.assertEqual(_read_gz_rows(self.outdir / "sample_snp_alleles.tsv.gz")[0][0], "sample_id")
        self.assertEqual(_read_gz_rows(self.outdir / "vcf_table.tsv.gz")[0], ["sample_id", "pos", "alt"])

    def test_sites_deduplicated_across_shards(self):
        rows = _read_gz_rows(self.outdir / "snp_sites.tsv.gz")[1:]
        self.assertEqual(len(rows), 8)
        keys = [row[0] for row in rows]
        self.assertEqual(len(keys), len(set(keys)))

    def test_sites_keep_max_qual_across_shards(self):
        rows = _read_gz_rows(self.outdir / "snp_sites.tsv.gz")[1:]
        by_pos = {int(row[2]): float(row[5]) for row in rows}
        # У SAMPLE2 QUAL втрое выше, значит должен победить он
        self.assertAlmostEqual(by_pos[1977], 1643.47 * 3, places=1)

    def test_alleles_and_vcf_rows_concatenated(self):
        self.assertEqual(len(_read_gz_rows(self.outdir / "sample_snp_alleles.tsv.gz")) - 1, 16)
        self.assertEqual(len(_read_gz_rows(self.outdir / "vcf_table.tsv.gz")) - 1, 14)

    def test_manifest(self):
        manifest = dict(
            line.split("\t")
            for line in (self.outdir / "snp_db_manifest.tsv").read_text().splitlines()
        )
        self.assertEqual(manifest["shards"], "2")
        self.assertEqual(manifest["samples"], "2")
        self.assertEqual(manifest["sites_rows_before_dedup"], "16")
        self.assertEqual(manifest["sites_rows"], "8")
        self.assertEqual(manifest["alleles_rows"], "16")
        self.assertEqual(manifest["vcf_rows"], "14")
```

- [ ] **Step 2: Запустить тест и убедиться, что он падает**

Run: `cd /home/zerg/git/tb-lite && python3 -m unittest tests.test_merge_snp_shards -v`
Expected: FAIL — `bash: .../bin/merge_snp_shards.sh: No such file or directory`

- [ ] **Step 3: Реализовать скрипт склейки**

`bin/merge_snp_shards.sh`:

```bash
#!/usr/bin/env bash
# Склейка шардов SNP_DB_SHARDS в итоговые таблицы для импорта в TB Platform.
set -euo pipefail

SHARD_DIR=""
OUT_DIR=""
THREADS=4
SORT_MEM="2G"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --shard-dir) SHARD_DIR="$2"; shift 2 ;;
        --out-dir)   OUT_DIR="$2";   shift 2 ;;
        --threads)   THREADS="$2";   shift 2 ;;
        --sort-mem)  SORT_MEM="$2";  shift 2 ;;
        *) echo "Неизвестный аргумент: $1" >&2; exit 2 ;;
    esac
done

if [[ -z "$SHARD_DIR" || -z "$OUT_DIR" ]]; then
    echo "Использование: merge_snp_shards.sh --shard-dir DIR --out-dir DIR [--threads N] [--sort-mem 2G]" >&2
    exit 2
fi

mkdir -p "$OUT_DIR"
SORT_TMP="$(mktemp -d "${OUT_DIR}/sort_tmp.XXXXXX")"
trap 'rm -rf "$SORT_TMP"' EXIT

if command -v pigz >/dev/null 2>&1; then
    compress() { pigz -p "$THREADS"; }
else
    compress() { gzip; }
fi

SITES_HEADER=$'site_key\tchrom\tpos\tref\talt\tqual\teffect\timpact\tgene\tgene_id\tproduct_name\tstrand\tfeature_id\tbiotype\trank\thgvs_c\thgvs_p\tcds_pos\taa_pos'
ALLELES_HEADER=$'sample_id\tsite_key\tallele\tdp\taf\tqual\tfilter'
VCF_HEADER=$'sample_id\tpos\talt'

shopt -s nullglob
site_shards=( "$SHARD_DIR"/sites.*.tsv.gz )
allele_shards=( "$SHARD_DIR"/alleles.*.tsv.gz )
vcf_shards=( "$SHARD_DIR"/vcfrows.*.tsv.gz )
count_files=( "$SHARD_DIR"/counts.*.tsv )
sample_files=( "$SHARD_DIR"/samples.*.txt )
shopt -u nullglob

if (( ${#site_shards[@]} == 0 )); then
    echo "ERROR: в $SHARD_DIR нет шардов sites.*.tsv.gz" >&2
    exit 1
fi

# snp_sites: одна строка на site_key, побеждает максимальный QUAL (как в проде).
# Поле 6 — qual; пустое значение считаем минимально возможным.
sites_rows=$(
    {
        printf '%s\n' "$SITES_HEADER"
        zcat "${site_shards[@]}" \
            | LC_ALL=C sort -t $'\t' -k1,1 -S "$SORT_MEM" --parallel="$THREADS" -T "$SORT_TMP" \
            | awk -F'\t' -v OFS='\t' '
                function q(v) { return (v == "" ? -1e308 : v + 0) }
                NR == 1 { key = $1; line = $0; best = q($6); next }
                $1 == key { if (q($6) > best) { line = $0; best = q($6) } ; next }
                { print line; key = $1; line = $0; best = q($6) }
                END { if (NR > 0) print line }
              '
    } | tee >(tail -n +2 | wc -l > "$SORT_TMP/sites_rows") | compress > "$OUT_DIR/snp_sites.tsv.gz"
    cat "$SORT_TMP/sites_rows"
)

{ printf '%s\n' "$ALLELES_HEADER"; zcat "${allele_shards[@]}"; } | compress > "$OUT_DIR/sample_snp_alleles.tsv.gz"
{ printf '%s\n' "$VCF_HEADER";     zcat "${vcf_shards[@]}";    } | compress > "$OUT_DIR/vcf_table.tsv.gz"

cat "${sample_files[@]}" | LC_ALL=C sort -u > "$OUT_DIR/samples.txt"

sum_metric() {
    local metric="$1"
    awk -F'\t' -v m="$metric" '$1 == m { total += $2 } END { print total + 0 }' "${count_files[@]}"
}

{
    printf 'shards\t%s\n'                  "${#site_shards[@]}"
    printf 'samples\t%s\n'                 "$(wc -l < "$OUT_DIR/samples.txt")"
    printf 'sites_rows_before_dedup\t%s\n' "$(sum_metric sites)"
    printf 'sites_rows\t%s\n'              "$sites_rows"
    printf 'alleles_rows\t%s\n'            "$(sum_metric alleles)"
    printf 'vcf_rows\t%s\n'                "$(sum_metric vcfrows)"
} > "$OUT_DIR/snp_db_manifest.tsv"

echo "OK: итоговые таблицы в $OUT_DIR"
```

Сделать исполняемым: `chmod +x bin/merge_snp_shards.sh`

- [ ] **Step 4: Запустить тест**

Run: `cd /home/zerg/git/tb-lite && python3 -m unittest tests.test_merge_snp_shards -v`
Expected: PASS, 6 тестов

- [ ] **Step 5: Коммит**

```bash
cd /home/zerg/git/tb-lite
git add bin/merge_snp_shards.sh tests/test_merge_snp_shards.py
git commit -m "feat(snp-db): склейка шардов в итоговые таблицы с дедупликацией по max(qual)"
```

---

### Task 4: Генератор `import_snp.sql`

**Files:**
- Create: `bin/write_import_snp_sql.py`
- Test: `tests/test_import_snp_sql.py` (только проверка текста; прогон против БД — Task 7)

**Interfaces:**
- Consumes: имена файлов из Task 3 (`snp_sites.tsv.gz`, `sample_snp_alleles.tsv.gz`, `vcf_table.tsv.gz`).
- Produces: CLI `write_import_snp_sql.py --data-dir <dir> -o <file> [--chunk-size N]`. Пишет самодостаточный psql-скрипт.

- [ ] **Step 1: Написать падающий тест**

`tests/test_import_snp_sql.py`:

```python
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
GENERATOR = REPO_ROOT / "bin" / "write_import_snp_sql.py"


class WriteImportSnpSqlTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        out = Path(cls._tmp.name) / "import_snp.sql"
        subprocess.run(
            [sys.executable, str(GENERATOR), "--data-dir", "/data/snp",
             "-o", str(out), "--chunk-size", "5000"],
            check=True,
        )
        cls.sql = out.read_text(encoding="utf-8")

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_stops_on_error(self):
        self.assertIn("\\set ON_ERROR_STOP on", self.sql)

    def test_copies_all_three_files(self):
        self.assertIn("/data/snp/snp_sites.tsv.gz", self.sql)
        self.assertIn("/data/snp/sample_snp_alleles.tsv.gz", self.sql)
        self.assertIn("/data/snp/vcf_table.tsv.gz", self.sql)

    def test_idempotent_inserts(self):
        self.assertIn("ON CONFLICT (site_key) DO NOTHING", self.sql)
        self.assertIn("ON CONFLICT (sample_id, site_id) DO NOTHING", self.sql)
        self.assertIn("ON CONFLICT (mask_key, sample_id) DO UPDATE", self.sql)

    def test_vcf_table_is_delete_then_insert(self):
        self.assertIn("DELETE FROM vcf_table v", self.sql)

    def test_only_samples_present_in_general(self):
        self.assertIn('FROM general g WHERE g."ID"', self.sql)

    def test_profiles_use_canonical_site_id(self):
        self.assertIn("min(site_id) OVER (PARTITION BY chrom, pos, ref, alt)", self.sql)

    def test_chunk_size_applied(self):
        self.assertIn("/ 5000", self.sql)

    def test_data_load_is_one_transaction_and_profiles_are_not(self):
        self.assertLess(self.sql.index("BEGIN;"), self.sql.index("COMMIT;"))
        self.assertGreater(self.sql.index("DO $do$"), self.sql.index("COMMIT;"))
```

- [ ] **Step 2: Запустить тест и убедиться, что он падает**

Run: `cd /home/zerg/git/tb-lite && python3 -m unittest tests.test_import_snp_sql -v`
Expected: FAIL — `can't open file .../bin/write_import_snp_sql.py`

- [ ] **Step 3: Реализовать генератор**

`bin/write_import_snp_sql.py`:

```python
#!/usr/bin/env python3
"""Генерация psql-скрипта импорта SNP-таблиц в TB Platform."""
from __future__ import annotations

import argparse
from pathlib import Path

TEMPLATE = r"""-- Импорт SNP-таблиц TB Platform. Сгенерировано пайплайном tb-lite.
--
-- Запуск (файлы должны быть видны тому процессу, который выполняет psql,
-- поэтому при запуске в контейнере каталог с данными нужно смонтировать):
--     psql "$DATABASE_URL" -f import_snp.sql
--
-- Порядок важен: general.tsv должен быть импортирован раньше, на
-- sample_snp_alleles.sample_id висит FK на general."ID".
-- Повторный запуск безопасен на каждом шаге.

\set ON_ERROR_STOP on
\timing on

BEGIN;

\echo '=== staging ==='

CREATE TEMP TABLE stg_sites (
    site_key text, chrom text, pos integer, ref text, alt text, qual double precision,
    effect text, impact text, gene text, gene_id text, product_name text, strand text,
    feature_id text, biotype text, rank text, hgvs_c text, hgvs_p text,
    cds_pos text, aa_pos text
);
\copy stg_sites FROM PROGRAM 'zcat ''{data_dir}/snp_sites.tsv.gz''' WITH (FORMAT csv, DELIMITER E'\t', HEADER true, NULL '', QUOTE E'\b')

CREATE TEMP TABLE stg_alleles (
    sample_id text, site_key text, allele text, dp integer, af double precision,
    qual double precision, filter text
);
\copy stg_alleles FROM PROGRAM 'zcat ''{data_dir}/sample_snp_alleles.tsv.gz''' WITH (FORMAT csv, DELIMITER E'\t', HEADER true, NULL '', QUOTE E'\b')

CREATE TEMP TABLE stg_vcf (sample_id text, pos integer, alt text);
\copy stg_vcf FROM PROGRAM 'zcat ''{data_dir}/vcf_table.tsv.gz''' WITH (FORMAT csv, DELIMITER E'\t', HEADER true, NULL '', QUOTE E'\b')

\echo '=== snp_sites ==='

INSERT INTO snp_sites (
    site_key, chrom, pos, ref, alt, qual, effect, impact, gene, gene_id,
    product_name, strand, feature_id, biotype, rank, hgvs_c, hgvs_p, cds_pos, aa_pos
)
SELECT
    site_key, chrom, pos, ref, alt, qual, effect, impact, gene, gene_id,
    product_name, strand, feature_id, biotype, rank, hgvs_c, hgvs_p, cds_pos, aa_pos
FROM stg_sites
ON CONFLICT (site_key) DO NOTHING;

CREATE INDEX ON stg_alleles (site_key);
ANALYZE stg_alleles;

\echo '=== sample_snp_alleles ==='

INSERT INTO sample_snp_alleles (sample_id, site_id, allele, dp, af, qual, filter)
SELECT a.sample_id, s.site_id, a.allele, a.dp, a.af, a.qual, a.filter
FROM stg_alleles a
JOIN snp_sites s ON s.site_key = a.site_key
WHERE EXISTS (SELECT 1 FROM general g WHERE g."ID" = a.sample_id)
ON CONFLICT (sample_id, site_id) DO NOTHING;

\echo '=== vcf_table ==='

CREATE TEMP TABLE stg_vcf_samples AS SELECT DISTINCT sample_id FROM stg_vcf;

DELETE FROM vcf_table v
USING stg_vcf_samples s
WHERE v.sample_id = s.sample_id;

INSERT INTO vcf_table (sample_id, pos, alt)
SELECT v.sample_id, v.pos, v.alt
FROM stg_vcf v
WHERE EXISTS (SELECT 1 FROM general g WHERE g."ID" = v.sample_id);

CREATE TEMP TABLE stg_new_samples AS
SELECT sample_id, ((row_number() OVER (ORDER BY sample_id)) - 1) / {chunk_size} AS chunk_no
FROM (
    SELECT DISTINCT a.sample_id
    FROM stg_alleles a
    WHERE EXISTS (SELECT 1 FROM general g WHERE g."ID" = a.sample_id)
) t;

CREATE INDEX ON stg_new_samples (chunk_no);
ANALYZE stg_new_samples;

COMMIT;

\echo '=== sample_snp_profiles ==='

-- Отдельно от загрузки данных: падение на профилях не должно откатывать
-- snp_sites/sample_snp_alleles/vcf_table. COMMIT внутри DO допустим, потому
-- что блок выполняется вне явной транзакции.
DO $do$
DECLARE
    masks text[];
    mask text;
    chunk integer;
    max_chunk integer;
BEGIN
    SELECT array_agg(mask_key ORDER BY sort_order) INTO masks FROM snp_masks;
    IF masks IS NULL THEN
        RAISE NOTICE 'Таблица snp_masks пуста — профили не строятся. Загрузите маски (deploy/rebuild_snp_profiles_masked.sql).';
        RETURN;
    END IF;

    SELECT max(chunk_no) INTO max_chunk FROM stg_new_samples;
    IF max_chunk IS NULL THEN
        RAISE NOTICE 'Нет новых образцов — профили не строятся.';
        RETURN;
    END IF;

    FOREACH mask IN ARRAY masks LOOP
        DROP TABLE IF EXISTS tmp_canon;
        CREATE TEMP TABLE tmp_canon AS
        WITH snp AS (
            SELECT ss.site_id, ss.chrom, ss.pos, ss.ref, ss.alt
            FROM snp_sites ss
            WHERE char_length(ss.ref) = 1 AND char_length(ss.alt) = 1
              AND upper(ss.ref) IN ('A','C','G','T','N')
              AND upper(ss.alt) IN ('A','C','G','T','N')
              AND (
                    mask = 'none'
                    OR NOT EXISTS (
                        SELECT 1 FROM snp_mask_regions r
                        WHERE r.mask_key = mask AND r.pos_range @> ss.pos
                    )
              )
        )
        SELECT site_id, min(site_id) OVER (PARTITION BY chrom, pos, ref, alt) AS cid
        FROM snp;

        CREATE INDEX ON tmp_canon (site_id);
        ANALYZE tmp_canon;

        FOR chunk IN 0..max_chunk LOOP
            INSERT INTO sample_snp_profiles (mask_key, sample_id, site_ids, snp_count)
            SELECT mask,
                   a.sample_id,
                   array_agg(DISTINCT c.cid ORDER BY c.cid),
                   count(DISTINCT c.cid)
            FROM sample_snp_alleles a
            JOIN tmp_canon c ON c.site_id = a.site_id
            WHERE a.sample_id IN (
                SELECT sample_id FROM stg_new_samples WHERE chunk_no = chunk
            )
            GROUP BY a.sample_id
            ON CONFLICT (mask_key, sample_id) DO UPDATE
                SET site_ids = EXCLUDED.site_ids,
                    snp_count = EXCLUDED.snp_count,
                    updated_at = now();
            COMMIT;
        END LOOP;
    END LOOP;

    DROP TABLE IF EXISTS tmp_canon;
END
$do$;

\echo '=== ANALYZE ==='
ANALYZE snp_sites;
ANALYZE sample_snp_profiles;

\echo 'Готово. vcf_table и sample_snp_alleles стоит проанализировать отдельно:'
\echo '  ANALYZE vcf_table; ANALYZE sample_snp_alleles;'
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True,
                        help="Каталог с snp_sites.tsv.gz / sample_snp_alleles.tsv.gz / vcf_table.tsv.gz "
                             "в той файловой системе, откуда будет запущен psql.")
    parser.add_argument("-o", "--out", required=True, type=Path)
    parser.add_argument("--chunk-size", type=int, default=5000,
                        help="Сколько образцов обрабатывать за одну транзакцию при построении профилей.")
    args = parser.parse_args()

    data_dir = args.data_dir.rstrip("/")
    args.out.write_text(
        TEMPLATE.format(data_dir=data_dir, chunk_size=args.chunk_size),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Запустить тест**

Run: `cd /home/zerg/git/tb-lite && python3 -m unittest tests.test_import_snp_sql -v`
Expected: PASS, 8 тестов

- [ ] **Step 5: Коммит**

```bash
cd /home/zerg/git/tb-lite
git add bin/write_import_snp_sql.py tests/test_import_snp_sql.py
git commit -m "feat(snp-db): генератор import_snp.sql со staging и инкрементальными профилями"
```

---

### Task 5: Nextflow-процесс `SNP_DB_SHARDS`

**Files:**
- Create: `modules/local/snp_db/shards/main.nf`
- Create: `modules/local/snp_db/shards/environment.yml`
- Create: `subworkflows/local/snp_db.nf`
- Create: `tests/nf/snp_db_shards_smoke.nf`
- Modify: `workflows/tblite.nf`
- Modify: `nextflow.config`
- Modify: `nextflow_schema.json`
- Modify: `conf/modules.config`

**Interfaces:**
- Consumes: `bin/vcf_to_snp_shards.py` из Task 2; канал `CALLVAR.out.ann` — `tuple(sample_id, ann_vcf)`.
- Produces: workflow `SNP_DB(ann_vcfs)`, где `ann_vcfs` — канал путей к `*.annotated.ann.vcf`. Публикует шарды в `${params.outdir}/snp_db/shards`.

- [ ] **Step 1: Добавить параметр `skip_snp_db`**

В `nextflow.config` в блок `// Reporting` после `skip_snp_matrix`:

```groovy
    skip_snp_db         = false
```

В `nextflow_schema.json` — рядом с описанием `skip_snp_matrix` добавить:

```json
"skip_snp_db": {
    "type": "boolean",
    "default": false,
    "description": "Не генерировать шарды SNP-таблиц (snp_sites / sample_snp_alleles / vcf_table) для TB Platform."
}
```

`params.batch_tag` уже существует и уже передаётся из `run_batches.sh:413` — менять его не нужно.

- [ ] **Step 2: Написать модуль процесса**

`modules/local/snp_db/shards/main.nf`:

```groovy
process SNP_DB_SHARDS {
    tag "${shard_tag}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tb-lite/tb-platform-tables:1.1' :
        'tb-lite/tb-platform-tables:1.1' }"

    input:
    tuple val(shard_tag), path(ann_vcfs)

    output:
    path "sites.${shard_tag}.tsv.gz",    emit: sites
    path "alleles.${shard_tag}.tsv.gz",  emit: alleles
    path "vcfrows.${shard_tag}.tsv.gz",  emit: vcfrows
    path "samples.${shard_tag}.txt",     emit: samples
    path "counts.${shard_tag}.tsv",      emit: counts

    script:
    """
    : > vcfs.list
    for vcf_path in ${ann_vcfs}; do
        printf '%s\\n' "\$(readlink -f "\$vcf_path")" >> vcfs.list
    done

    python ${projectDir}/bin/vcf_to_snp_shards.py \\
        --file-list vcfs.list \\
        --tag ${shard_tag} \\
        -o . \\
        -j ${task.cpus}
    """
}
```

`modules/local/snp_db/shards/environment.yml`:

```yaml
---
channels:
  - conda-forge
dependencies:
  - conda-forge::python=3.12
```

- [ ] **Step 3: Написать subworkflow**

`subworkflows/local/snp_db.nf`:

```groovy
include { SNP_DB_SHARDS } from '../../modules/local/snp_db/shards/main'

// Детерминированный тег по составу пачки: перезапуск перезаписывает тот же
// файл, а не создаёт второй с дублями. Имена образцов в пачках не пересекаются,
// поэтому первого имени и размера достаточно для уникальности.
def snpDbChunkTag(List paths, String prefix) {
    def names = paths.collect { it.name }.sort()
    // tokenize('.')[0] вместо регулярки: в slashy-строке Groovy символ '$'
    // перед '/' даёт ошибку компиляции.
    return "${prefix}_${names.first().tokenize('.')[0]}_${names.size()}"
}

workflow SNP_DB {
    take:
        ann_vcfs   // канал путей к *.annotated.ann.vcf

    main:
        // run_batches.sh уже передаёт --batch_tag batch_${i}; без него тег
        // выводится из состава батча.
        tagged = ann_vcfs
            .collect()
            .map { paths ->
                tuple(params.batch_tag ? params.batch_tag as String
                                       : snpDbChunkTag(paths, 'batch'), paths)
            }

        SNP_DB_SHARDS(tagged)

    emit:
        sites   = SNP_DB_SHARDS.out.sites
        alleles = SNP_DB_SHARDS.out.alleles
        vcfrows = SNP_DB_SHARDS.out.vcfrows
}
```

Процесс принимает один tuple-канал, поэтому тег и список файлов не могут
разъехаться по порядку — в отличие от двух независимых входных каналов.

- [ ] **Step 4: Настроить publishDir**

В `conf/modules.config` внутри блока `process { ... }` добавить:

```groovy
    withName: 'SNP_DB_SHARDS' {
        publishDir = [
            path: { "${params.outdir}/snp_db/shards" },
            mode: params.mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
```

- [ ] **Step 5: Подключить в основной workflow**

В `workflows/tblite.nf` добавить импорт после строки с `VCF_ANNOTATION`:

```groovy
include { SNP_DB }      from '../subworkflows/local/snp_db'
```

И сразу после блока `if (!params.skip_snp_matrix) { ANN_TABLE(...) }`:

```groovy
        if (!params.skip_snp_db) {
            SNP_DB(callvar.ann.map { _sample_name, ann_vcf -> ann_vcf })
        }
```

- [ ] **Step 6: Написать smoke-workflow и прогнать процесс**

`tests/nf/snp_db_shards_smoke.nf`:

```groovy
nextflow.enable.dsl = 2

include { SNP_DB } from '../../subworkflows/local/snp_db'

workflow {
    SNP_DB(Channel.fromPath(params.smoke_vcfs, checkIfExists: true))
}
```

Run:

```bash
cd /home/zerg/git/tb-lite
rm -rf /tmp/snp_db_smoke && mkdir -p /tmp/snp_db_smoke
nextflow run tests/nf/snp_db_shards_smoke.nf \
    -profile docker \
    -w /tmp/snp_db_smoke/work \
    --outdir /tmp/snp_db_smoke/results \
    --batch_tag batch_smoke \
    --smoke_vcfs '/mnt/data1/tb-lite-runs/bc1286ed-5e2e-42db-a06c-0a5a02258824/results/annotate_vcf/*/*.annotated.ann.vcf'
cat /tmp/snp_db_smoke/results/snp_db/shards/counts.batch_smoke.tsv
```

Expected: процесс завершается успешно, `counts.batch_smoke.tsv` содержит `samples 3`, а `sites`/`alleles`/`vcfrows` — ненулевые значения.

- [ ] **Step 7: Коммит**

```bash
cd /home/zerg/git/tb-lite
git add modules/local/snp_db/shards/ subworkflows/local/snp_db.nf tests/nf/snp_db_shards_smoke.nf \
        workflows/tblite.nf nextflow.config nextflow_schema.json conf/modules.config
git commit -m "feat(snp-db): процесс SNP_DB_SHARDS и подключение к основному workflow"
```

---

### Task 6: Nextflow-процесс `SNP_DB_TABLES` в `batch_reports.nf`

**Files:**
- Create: `modules/local/snp_db/tables/main.nf`
- Create: `modules/local/snp_db/tables/environment.yml`
- Create: `tests/nf/snp_db_final_smoke.nf`
- Modify: `subworkflows/local/snp_db.nf`
- Modify: `batch_reports.nf`
- Modify: `conf/modules.config`

**Interfaces:**
- Consumes: `bin/merge_snp_shards.sh` (Task 3), `bin/write_import_snp_sql.py` (Task 4), `SNP_DB_SHARDS` (Task 5).
- Produces: workflow `SNP_DB_FINAL()` — сам находит шарды или строит их из `annotate_vcf/`; публикует итог в `${params.outdir}/Reports/tb-platform/snp`.

- [ ] **Step 1: Написать модуль `SNP_DB_TABLES`**

`modules/local/snp_db/tables/main.nf`:

```groovy
process SNP_DB_TABLES {
    tag "SNP DB tables"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tb-lite/tb-platform-tables:1.1' :
        'tb-lite/tb-platform-tables:1.1' }"

    input:
    path shards

    output:
    path "snp_sites.tsv.gz",           emit: sites
    path "sample_snp_alleles.tsv.gz",  emit: alleles
    path "vcf_table.tsv.gz",           emit: vcf_table
    path "samples.txt",                emit: samples
    path "snp_db_manifest.tsv",        emit: manifest
    path "import_snp.sql",             emit: import_sql

    script:
    // Шарды уже разложены Nextflow'ом в рабочем каталоге задачи, поэтому
    // --shard-dir . Итоговые имена (snp_sites.*, samples.txt) под маски
    // шардов (sites.*.tsv.gz, samples.*.txt) не попадают.
    def publish_dir = file(params.outdir).toAbsolutePath().toString() + '/Reports/tb-platform/snp'
    """
    bash ${projectDir}/bin/merge_snp_shards.sh \\
        --shard-dir . \\
        --out-dir . \\
        --threads ${task.cpus}

    # Путь в SQL — это место публикации, откуда файлы будет читать psql.
    python ${projectDir}/bin/write_import_snp_sql.py \\
        --data-dir "${publish_dir}" \\
        -o import_snp.sql
    """
}
```

`modules/local/snp_db/tables/environment.yml`:

```yaml
---
channels:
  - conda-forge
dependencies:
  - conda-forge::python=3.12
  - conda-forge::pigz
```

- [ ] **Step 2: Добавить `SNP_DB_FINAL` в subworkflow**

Дописать в `subworkflows/local/snp_db.nf`:

Добавить в начало файла ещё два импорта:

```groovy
include { SNP_DB_TABLES } from '../../modules/local/snp_db/tables/main'
include { SNP_DB_SHARDS as REBUILD_SHARDS } from '../../modules/local/snp_db/shards/main'
```

И сам workflow:

```groovy
workflow SNP_DB_FINAL {
    main:
        // Наличие шардов — это состояние файловой системы на момент запуска,
        // а не поток данных. Решаем в Groovy: так каждый канал потребляется
        // ровно один раз и ветки не конфликтуют.
        def shard_dir = file("${params.outdir}/snp_db/shards")
        def existing_shards = shard_dir.exists() ? files("${shard_dir}/sites.*.tsv.gz") : []

        if (existing_shards) {
            log.info("SNP DB: найдено ${existing_shards.size()} шардов в ${shard_dir}")
            shards = Channel
                .fromPath("${shard_dir}/{sites,alleles,vcfrows,samples,counts}.*")
                .collect()
        }
        else {
            log.info("SNP DB: шардов нет, собираю их из ${params.outdir}/annotate_vcf")
            chunks = Channel
                .fromPath("${params.outdir}/annotate_vcf/*/*.annotated.ann.vcf", checkIfExists: true)
                .collate(500)
                .map { chunk -> tuple(snpDbChunkTag(chunk, 'rebuilt'), chunk) }

            REBUILD_SHARDS(chunks)

            shards = REBUILD_SHARDS.out.sites
                .mix(
                    REBUILD_SHARDS.out.alleles,
                    REBUILD_SHARDS.out.vcfrows,
                    REBUILD_SHARDS.out.samples,
                    REBUILD_SHARDS.out.counts
                )
                .collect()
        }

        SNP_DB_TABLES(shards)

    emit:
        manifest = SNP_DB_TABLES.out.manifest
}
```

В ветке пересборки `params.batch_tag` намеренно игнорируется: пачек несколько,
и общий тег на все привёл бы к тому, что они затирают файлы друг друга.

- [ ] **Step 3: Подключить в `batch_reports.nf`**

В `batch_reports.nf` добавить импорт после строки с `REPORTS`:

```groovy
include { SNP_DB_FINAL } from './subworkflows/local/snp_db'
```

И в конце `workflow BATCH_REPORTS` после вызова `REPORTS(...)`:

```groovy
        if (!params.skip_snp_db) {
            SNP_DB_FINAL()
        }
```

- [ ] **Step 4: Настроить publishDir**

В `conf/modules.config` добавить:

```groovy
    withName: 'SNP_DB_TABLES' {
        publishDir = [
            path: { "${params.outdir}/Reports/tb-platform/snp" },
            mode: params.mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
    withName: 'REBUILD_SHARDS' {
        publishDir = [
            path: { "${params.outdir}/snp_db/shards" },
            mode: params.mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
```

- [ ] **Step 5: Прогнать ветку пересборки на реальном каталоге без шардов**

Отдельный smoke-workflow нужен потому, что `BATCH_REPORTS` объявляет входные каналы с `checkIfExists: true` и упадёт раньше `SNP_DB_FINAL`, если в каталоге прогона пуст хотя бы один из его глобов. Полный прогон `batch_reports.nf` проверяется в Task 8.

`tests/nf/snp_db_final_smoke.nf`:

```groovy
nextflow.enable.dsl = 2

include { SNP_DB_FINAL } from '../../subworkflows/local/snp_db'

workflow {
    SNP_DB_FINAL()
}
```

Run:

```bash
cd /home/zerg/git/tb-lite
rm -rf /tmp/snp_db_final && mkdir -p /tmp/snp_db_final
cp -r /mnt/data1/tb-lite-runs/bc1286ed-5e2e-42db-a06c-0a5a02258824/results /tmp/snp_db_final/results
rm -rf /tmp/snp_db_final/results/snp_db
nextflow run tests/nf/snp_db_final_smoke.nf \
    -profile docker \
    -w /tmp/snp_db_final/work \
    --outdir /tmp/snp_db_final/results
ls -la /tmp/snp_db_final/results/Reports/tb-platform/snp/
cat /tmp/snp_db_final/results/Reports/tb-platform/snp/snp_db_manifest.tsv
zcat /tmp/snp_db_final/results/Reports/tb-platform/snp/snp_sites.tsv.gz | head -2
grep data-dir -n /tmp/snp_db_final/results/Reports/tb-platform/snp/import_snp.sql || \
  grep -m1 'snp_sites.tsv.gz' /tmp/snp_db_final/results/Reports/tb-platform/snp/import_snp.sql
```

Expected: в `Reports/tb-platform/snp/` шесть файлов; в манифесте `shards 1` и `samples 3`; первая строка `snp_sites.tsv.gz` — заголовок `site_key ...`; в `import_snp.sql` путь абсолютный и указывает на `/tmp/snp_db_final/results/Reports/tb-platform/snp`.

- [ ] **Step 6: Коммит**

```bash
cd /home/zerg/git/tb-lite
git add modules/local/snp_db/tables/ subworkflows/local/snp_db.nf batch_reports.nf \
        conf/modules.config tests/nf/snp_db_final_smoke.nf
git commit -m "feat(snp-db): финальная сборка SNP-таблиц в batch_reports"
```

---

### Task 7: Интеграционный тест импорта в PostgreSQL

**Files:**
- Modify: `tests/test_import_snp_sql.py`

**Interfaces:**
- Consumes: артефакты Task 6 и `bin/write_import_snp_sql.py`.
- Produces: тест-класс `ImportSnpSqlIntegrationTest`, который прогоняет `import_snp.sql` против контейнера `tb_platform_db` и проверяет счётчики и идемпотентность.

- [ ] **Step 1: Написать падающий тест**

Дописать в `tests/test_import_snp_sql.py`:

```python
import shutil

DB_CONTAINER = "tb_platform_db"
PSQL = ["docker", "exec", "-e", "PGPASSWORD=tb_password", "-i", DB_CONTAINER,
        "psql", "-U", "tb_user", "-d", "tb_database", "-v", "ON_ERROR_STOP=1"]


def _docker_available():
    if shutil.which("docker") is None:
        return False
    probe = subprocess.run(["docker", "exec", DB_CONTAINER, "true"], capture_output=True)
    return probe.returncode == 0


def _sql(statement):
    done = subprocess.run(PSQL + ["-tAc", statement], capture_output=True, text=True)
    if done.returncode != 0:
        raise AssertionError(done.stderr)
    return done.stdout.strip()


@unittest.skipUnless(_docker_available(), "нужен запущенный контейнер tb_platform_db")
class ImportSnpSqlIntegrationTest(unittest.TestCase):
    SAMPLES = ("ERR4797591", "SAMPLE2")

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        tmp = Path(cls._tmp.name)

        # Данные готовим тем же путём, что и пайплайн: конвертер -> склейка.
        second = tmp / "SAMPLE2.annotated.ann.vcf"
        fixture = REPO_ROOT / "tests" / "data" / "ERR4797591.mini.ann.vcf"
        lines = []
        for line in fixture.read_text(encoding="utf-8").splitlines():
            lines.append(line.replace("ERR4797591", "SAMPLE2") if line.startswith("#CHROM") else line)
        second.write_text("\n".join(lines) + "\n", encoding="utf-8")

        shard_dir = tmp / "shards"
        for tag, vcf in (("batch_1", fixture), ("batch_2", second)):
            listing = tmp / f"{tag}.txt"
            listing.write_text(f"{vcf}\n", encoding="utf-8")
            subprocess.run(
                [sys.executable, str(REPO_ROOT / "bin" / "vcf_to_snp_shards.py"),
                 "--file-list", str(listing), "--tag", tag, "-o", str(shard_dir), "-j", "1"],
                check=True,
            )

        cls.data_dir = tmp / "data"
        subprocess.run(
            ["bash", str(REPO_ROOT / "bin" / "merge_snp_shards.sh"),
             "--shard-dir", str(shard_dir), "--out-dir", str(cls.data_dir), "--threads", "1"],
            check=True,
        )

        cls.container_dir = "/tmp/snp_import_test"
        subprocess.run(["docker", "exec", DB_CONTAINER, "rm", "-rf", cls.container_dir], check=True)
        subprocess.run(["docker", "cp", str(cls.data_dir), f"{DB_CONTAINER}:{cls.container_dir}"], check=True)

        subprocess.run(
            [sys.executable, str(REPO_ROOT / "bin" / "write_import_snp_sql.py"),
             "--data-dir", cls.container_dir, "-o", str(tmp / "import_snp.sql")],
            check=True,
        )
        cls.sql_path = tmp / "import_snp.sql"

        cls._seed_db()

    @classmethod
    def _seed_db(cls):
        ids = ", ".join(f"('{name}')" for name in cls.SAMPLES)
        _sql(f'INSERT INTO general ("ID") VALUES {ids} ON CONFLICT DO NOTHING;')
        _sql("""
            INSERT INTO snp_masks (mask_key, name, sort_order) VALUES
              ('none', 'Без маски', 0),
              ('rlc_lowmap', 'RLC + low-map', 1),
              ('dr_genes', 'DR-гены', 2),
              ('union', 'Объединение', 3)
            ON CONFLICT DO NOTHING;
        """)
        # Одна регион-маска, накрывающая pos 4013: профиль union должен стать короче none.
        _sql("""
            INSERT INTO snp_mask_regions (mask_key, chrom, pos_range) VALUES
              ('union', 'NC_000962.3', int4range(4000, 4100))
            ON CONFLICT DO NOTHING;
        """)

    @classmethod
    def tearDownClass(cls):
        ids = ", ".join(f"'{name}'" for name in cls.SAMPLES)
        _sql(f"DELETE FROM sample_snp_profiles WHERE sample_id IN ({ids});")
        _sql(f"DELETE FROM vcf_table WHERE sample_id IN ({ids});")
        _sql(f'DELETE FROM general WHERE "ID" IN ({ids});')
        _sql("DELETE FROM sample_snp_alleles;")
        _sql("DELETE FROM snp_sites;")
        _sql("DELETE FROM snp_mask_regions;")
        _sql("DELETE FROM snp_masks;")
        subprocess.run(["docker", "exec", DB_CONTAINER, "rm", "-rf", cls.container_dir], check=True)
        cls._tmp.cleanup()

    def _run_import(self):
        with open(self.sql_path, "rb") as handle:
            done = subprocess.run(PSQL + ["-f", "-"], stdin=handle, capture_output=True, text=True)
        self.assertEqual(done.returncode, 0, done.stderr)
        return done.stdout

    def _counts(self):
        return {
            "snp_sites": int(_sql("SELECT count(*) FROM snp_sites;")),
            "sample_snp_alleles": int(_sql("SELECT count(*) FROM sample_snp_alleles;")),
            "vcf_table": int(_sql("SELECT count(*) FROM vcf_table;")),
            "sample_snp_profiles": int(_sql("SELECT count(*) FROM sample_snp_profiles;")),
        }

    def test_import_then_reimport_is_idempotent(self):
        self._run_import()
        first = self._counts()

        manifest = dict(
            line.split("\t")
            for line in (self.data_dir / "snp_db_manifest.tsv").read_text().splitlines()
        )
        self.assertEqual(first["snp_sites"], int(manifest["sites_rows"]))
        self.assertEqual(first["sample_snp_alleles"], int(manifest["alleles_rows"]))
        self.assertEqual(first["vcf_table"], int(manifest["vcf_rows"]))
        # 2 образца x 4 маски
        self.assertEqual(first["sample_snp_profiles"], 8)

        self._run_import()
        self.assertEqual(self._counts(), first)

    def test_mask_shrinks_profile(self):
        self._run_import()
        none_count = int(_sql(
            "SELECT snp_count FROM sample_snp_profiles "
            "WHERE mask_key='none' AND sample_id='ERR4797591';"
        ))
        union_count = int(_sql(
            "SELECT snp_count FROM sample_snp_profiles "
            "WHERE mask_key='union' AND sample_id='ERR4797591';"
        ))
        self.assertGreater(none_count, union_count)
```

- [ ] **Step 2: Запустить тест и убедиться, что он падает**

Run: `cd /home/zerg/git/tb-lite && python3 -m unittest tests.test_import_snp_sql.ImportSnpSqlIntegrationTest -v`
Expected: FAIL — падение на `_run_import` или на несовпадении счётчиков (если в `import_snp.sql` есть ошибка синтаксиса или логики).

Если тест **проходит сразу** — это тоже валидный результат: значит SQL из Task 4 корректен. Тогда переходить к Step 4.

- [ ] **Step 3: Починить `bin/write_import_snp_sql.py` по тексту ошибки**

Типичные места: `COMMIT` внутри `DO` требует, чтобы блок выполнялся вне явной транзакции (после `COMMIT;`); `FOREACH ... IN ARRAY` не должен быть заменён на `FOR ... IN SELECT`, иначе `COMMIT` упадёт с «cannot commit while a portal is pinned»; имя переменной `mask` не должно конфликтовать с колонкой.

- [ ] **Step 4: Запустить весь набор тестов**

Run: `cd /home/zerg/git/tb-lite && python3 -m unittest discover -s tests -t . -v`
Expected: PASS, все тесты

- [ ] **Step 5: Коммит**

```bash
cd /home/zerg/git/tb-lite
git add tests/test_import_snp_sql.py bin/write_import_snp_sql.py
git commit -m "test(snp-db): интеграционный тест импорта и идемпотентности"
```

---

### Task 8: Сквозная проверка на реальном прогоне и документация

**Files:**
- Modify: `README.md`
- Modify: `CHANGELOG.md`

**Interfaces:**
- Consumes: всё предыдущее.
- Produces: раздел README с инструкцией по импорту.

- [ ] **Step 1: Прогнать пайплайн целиком на трёх образцах**

Это единственная проверка, которая проходит по тому же коду, что и продовый запуск: `main.nf` → `CALLVAR` → `SNP_DB`.

```bash
cd /home/zerg/git/tb-lite
rm -rf /tmp/snp_db_e2e && mkdir -p /tmp/snp_db_e2e
nextflow run main.nf \
    -profile docker \
    -w /tmp/snp_db_e2e/work \
    --input /mnt/data1/tb-lite-runs/bc1286ed-5e2e-42db-a06c-0a5a02258824/samplesheet.host.csv \
    --outdir /tmp/snp_db_e2e/results \
    --batch_tag batch_1 \
    --skip_kraken --skip_snp_matrix
ls /tmp/snp_db_e2e/results/snp_db/shards/
cat /tmp/snp_db_e2e/results/snp_db/shards/counts.batch_1.tsv
```

Expected: шарды с тегом `batch_1`, `samples 3`.

- [ ] **Step 2: Собрать итог тем же прогоном batch_reports**

```bash
cd /home/zerg/git/tb-lite
nextflow run batch_reports.nf \
    -profile docker \
    -w /tmp/snp_db_e2e/work_reports \
    --outdir /tmp/snp_db_e2e/results \
    --skip_multiqc
cat /tmp/snp_db_e2e/results/Reports/tb-platform/snp/snp_db_manifest.tsv
```

Expected: манифест с `shards 1`, `samples 3`; шарды переиспользованы, а не пересобраны из `annotate_vcf/`.

- [ ] **Step 3: Сверить site_key с продом на реальных данных**

```bash
zcat /tmp/snp_db_e2e/results/Reports/tb-platform/snp/snp_sites.tsv.gz \
  | awk -F'\t' '$3==1977 || $3==4013 {print $3"\t"$1}'
```

Expected:

```
1977	7229c7508aa7fb285706a91c398ca4b22bdccee9
4013	f03d041b126f9db3c71d4f62fb698395201f8ee5
```

Если значения разошлись — реализация несовместима с существующей базой, дальше идти нельзя.

- [ ] **Step 4: Дописать README**

Добавить в `README.md` раздел:

```markdown
## SNP-таблицы для TB Platform

Пайплайн генерирует файлы для добавления образцов на сайт. Итог — в
`<outdir>/Reports/tb-platform/snp/`:

| Файл | Назначение |
|---|---|
| `snp_sites.tsv.gz` | каталог аннотированных сайтов (`site_id` выдаёт БД по `site_key`) |
| `sample_snp_alleles.tsv.gz` | аллели образцов; нужны для SNP matrix TSV и поиска по SNP |
| `vcf_table.tsv.gz` | плоская таблица `sample_id/pos/alt` для поиска и clade-признаков |
| `import_snp.sql` | psql-скрипт загрузки |
| `snp_db_manifest.tsv` | счётчики строк для сверки |
| `samples.txt` | список загружаемых образцов |

Пишется всё содержимое VCF: SNP, инделы, complex и mnp, без фильтров и масок.

Отключается флагом `--skip_snp_db`.

### Импорт

`general.tsv` должен быть импортирован раньше — на `sample_snp_alleles.sample_id`
висит внешний ключ на `general."ID"`.

```bash
psql "$DATABASE_URL" -f <outdir>/Reports/tb-platform/snp/import_snp.sql
```

`\copy` выполняется на стороне клиента, поэтому каталог с `.tsv.gz` должен быть
доступен тому процессу, который запускает psql. При запуске psql в контейнере
каталог нужно смонтировать, а `import_snp.sql` перегенерировать с путём внутри
контейнера:

```bash
python bin/write_import_snp_sql.py --data-dir /mnt/snp -o import_snp.sql
```

Повторный запуск импорта безопасен: `snp_sites` схлопывается по `site_key`,
`sample_snp_alleles` — по первичному ключу, `vcf_table` перезаписывается по
списку образцов, профили обновляются.

Профили `sample_snp_profiles` строятся последним шагом того же SQL по маскам из
`snp_masks`/`snp_mask_regions`. Если таблица масок пуста, шаг пропускается с
предупреждением — маски заливаются отдельно
(`tb-platform/deploy/rebuild_snp_profiles_masked.sql`).
```

- [ ] **Step 5: Дописать CHANGELOG**

Добавить в `CHANGELOG.md` в начало списка изменений:

```markdown
- Добавлена генерация SNP-таблиц для TB Platform: `snp_sites.tsv.gz`,
  `sample_snp_alleles.tsv.gz`, `vcf_table.tsv.gz` и `import_snp.sql`
  в `Reports/tb-platform/snp/`. Отключается флагом `--skip_snp_db`.
```

- [ ] **Step 6: Финальный прогон тестов и коммит**

Run: `cd /home/zerg/git/tb-lite && python3 -m unittest discover -s tests -t . -v`
Expected: PASS

```bash
cd /home/zerg/git/tb-lite
git add README.md CHANGELOG.md
git commit -m "docs(snp-db): описание SNP-таблиц и порядка импорта"
```

---

## Порядок и зависимости

```
Task 1 (ядро) ─> Task 2 (CLI/шарды) ─> Task 3 (склейка) ─┐
                                                          ├─> Task 6 (batch_reports) ─> Task 8 (e2e + docs)
                        Task 4 (SQL) ────────────────────┘
                                 │
                        Task 5 (SNP_DB_SHARDS) ──────────┘
                                 │
                        Task 7 (интеграция с БД)
```

Task 4 и Task 5 независимы друг от друга и могут идти в любом порядке после Task 3.
Task 7 требует Task 3 и Task 4, но не требует Nextflow-частей.
