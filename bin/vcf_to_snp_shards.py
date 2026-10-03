#!/usr/bin/env python3
"""Конвертация per-sample аннотированных VCF в шарды для таблиц TB Platform.

Пишет три шарда без заголовка:
  sites.<tag>.tsv.gz     site_key + атрибуты снипа (site_id раздаёт Postgres)
  alleles.<tag>.tsv.gz   образец x site_key
  vcfrows.<tag>.tsv.gz   образец x позиция x альтернативная аллель

Фильтрации нет: берём каждую ALT-аллель каждой записи VCF.

Логика извлечения перенесена без изменений из
tb-platform/new_temp_data/vcf_to_snp_db_tsv_parallel_scan.py: формула site_key
должна побитово совпадать с продовой таблицей snp_sites, иначе новые данные не
схлопнутся с существующими через ON CONFLICT (site_key).
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
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

SITE_COLUMNS = [
    "site_key", "chrom", "pos", "ref", "alt", "qual", "effect", "impact", "gene",
    "gene_id", "product_name", "strand", "feature_id", "biotype", "rank",
    "hgvs_c", "hgvs_p", "cds_pos", "aa_pos",
]
ALLELE_COLUMNS = ["sample_id", "site_key", "allele", "dp", "af", "qual", "filter"]
VCF_COLUMNS = ["sample_id", "pos", "alt"]

_QUAL_IDX = SITE_COLUMNS.index("qual")


ANN_SPLIT_RE = re.compile(r",(?=[^,]*\|)")
MISSING = {"", ".", "NA", "NaN", "nan", "null", "None"}


def eprint(*args: object) -> None:
    print(*args, file=sys.stderr)


def clean(value: object) -> str:
    if value is None:
        return ""
    s = str(value)
    if s in MISSING:
        return ""
    return s.replace("\t", " ").replace("\n", " ").replace("\r", " ")


def as_float(value: object) -> Optional[float]:
    s = clean(value)
    if not s:
        return None
    try:
        return float(s)
    except Exception:
        return None


def as_int(value: object) -> Optional[int]:
    s = clean(value)
    if not s:
        return None
    try:
        return int(float(s))
    except Exception:
        return None


def float_to_tsv(value: Optional[float]) -> str:
    if value is None:
        return ""
    if value != value:  # NaN
        return ""
    return format(value, ".10g")


def int_to_tsv(value: Optional[int]) -> str:
    return "" if value is None else str(value)


def open_text(path: Path):
    name = path.name.lower()
    if name.endswith(".gz") or name.endswith(".bgz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def open_out(path: Path, gzip_output: bool = False):
    if gzip_output:
        return gzip.open(path, "wt", encoding="utf-8", compresslevel=6, newline="")
    return open(path, "w", encoding="utf-8", newline="")


def strip_vcf_suffix(path: Path) -> str:
    name = path.name
    for suffix in [".vcf.gz", ".vcf.bgz", ".vcf", ".gz", ".bgz"]:
        if name.endswith(suffix):
            sample = name[: -len(suffix)]
            break
    else:
        sample = path.stem

    # Many TB-Lite outputs may look like SAMPLE__SAMPLE.vcf.gz.
    # Collapse this to SAMPLE so sample_id can match the main platform table.
    parts = sample.split("__")
    if len(parts) == 2 and parts[0] == parts[1]:
        return parts[0]
    return sample


def parse_info(info: str) -> Dict[str, str]:
    result: Dict[str, str] = {}
    if not info or info == ".":
        return result
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
        else:
            result[item] = "true"
    return result


def split_list_value(value: str) -> List[str]:
    if not value or value == ".":
        return []
    return value.split(",")


def select_list_value(value: str, index: int) -> str:
    parts = split_list_value(value)
    if not parts:
        return ""
    if 0 <= index < len(parts):
        return parts[index]
    return parts[0]


def parse_format(format_text: str, sample_text: str) -> Dict[str, str]:
    if not format_text or format_text == "." or not sample_text or sample_text == ".":
        return {}
    keys = format_text.split(":")
    vals = sample_text.split(":")
    if len(vals) < len(keys):
        vals += [""] * (len(keys) - len(vals))
    return dict(zip(keys, vals))


def gt_alt_indices(gt: str, alt_count: int) -> Optional[List[int]]:
    """Return 0-based ALT indices from GT. None means GT is absent/unusable."""
    gt = clean(gt)
    if not gt:
        return None
    gt = gt.split(":", 1)[0]
    if gt in MISSING:
        return []
    tokens = re.split(r"[|/]", gt)
    indices: List[int] = []
    for tok in tokens:
        tok = tok.strip()
        if not tok or tok == ".":
            continue
        try:
            allele_index = int(tok)
        except Exception:
            continue
        if allele_index <= 0:
            continue
        alt_idx = allele_index - 1
        if 0 <= alt_idx < alt_count:
            indices.append(alt_idx)
    return sorted(set(indices))


def estimate_af_from_depths(fmt: Dict[str, str], info: Dict[str, str], alt_idx: int, dp: Optional[int]) -> Optional[float]:
    # FreeBayes often uses AO/RO; GATK often uses AD.
    ao = clean(fmt.get("AO") or info.get("AO"))
    if ao:
        ao_value = as_float(select_list_value(ao, alt_idx))
        if ao_value is not None and dp and dp > 0:
            return ao_value / dp

    ad = clean(fmt.get("AD") or info.get("AD"))
    if ad:
        parts = [as_float(x) for x in ad.split(",")]
        if len(parts) >= alt_idx + 2:
            alt_depth = parts[alt_idx + 1]
            total = sum(x for x in parts if x is not None)
            if alt_depth is not None and total > 0:
                return alt_depth / total

    return None


def extract_metrics(fmt: Dict[str, str], info: Dict[str, str], alt_idx: int) -> Tuple[Optional[int], Optional[float]]:
    dp = as_int(fmt.get("DP") or info.get("DP"))

    af_text = clean(
        fmt.get("AF")
        or fmt.get("VAF")
        or fmt.get("VF")
        or info.get("AF")
        or info.get("VAF")
    )
    af = as_float(select_list_value(af_text, alt_idx)) if af_text else None

    if af is None:
        af = estimate_af_from_depths(fmt, info, alt_idx, dp)

    return dp, af


def ann_entries_for_alt(info: Dict[str, str], alt: str, single_alt: bool) -> List[List[str]]:
    ann_value = info.get("ANN", "")
    if not ann_value:
        return [[]]

    entries = ANN_SPLIT_RE.split(ann_value)
    parsed = [entry.split("|") for entry in entries if entry]

    selected: List[List[str]] = []
    for parts in parsed:
        allele = parts[0] if parts else ""
        if allele == alt:
            selected.append(parts)

    # In some indel/complex cases SnpEff's ANN allele may not be byte-identical to ALT.
    # If the record is single-ALT, using all ANN entries is usually less wrong than losing annotation.
    if not selected and single_alt:
        selected = parsed

    if not selected:
        return [[]]
    return selected


def ann_field(parts: Sequence[str], index: int) -> str:
    return clean(parts[index]) if len(parts) > index else ""


def site_row_from_ann(
    chrom: str,
    pos: str,
    ref: str,
    alt: str,
    qual: Optional[float],
    info: Dict[str, str],
    ann_parts: Sequence[str],
) -> List[str]:
    # SnpEff ANN layout:
    # 0 Allele
    # 1 Annotation / EFFECT
    # 2 Annotation_Impact
    # 3 Gene_Name
    # 4 Gene_ID
    # 5 Feature_Type
    # 6 Feature_ID
    # 7 Transcript_BioType
    # 8 Rank/Total
    # 9 HGVS.c
    # 10 HGVS.p
    # 11 cDNA.pos/cDNA.length
    # 12 CDS.pos/CDS.length
    # 13 AA.pos/AA.length
    product_name = clean(
        info.get("name")
        or info.get("Name")
        or info.get("product")
        or info.get("PRODUCT")
        or info.get("gene_name")
    )
    strand = clean(info.get("strand") or info.get("STRAND"))

    effect = ann_field(ann_parts, 1)
    impact = ann_field(ann_parts, 2)
    gene = ann_field(ann_parts, 3)
    gene_id = ann_field(ann_parts, 4)
    feature_id = ann_field(ann_parts, 6)
    biotype = ann_field(ann_parts, 7)
    rank = ann_field(ann_parts, 8)
    hgvs_c = ann_field(ann_parts, 9)
    hgvs_p = ann_field(ann_parts, 10)
    cds_pos = ann_field(ann_parts, 12)
    aa_pos = ann_field(ann_parts, 13)

    fields_without_key = [
        clean(chrom),
        clean(pos),
        clean(ref),
        clean(alt),
        float_to_tsv(qual),
        effect,
        impact,
        gene,
        gene_id,
        product_name,
        strand,
        feature_id,
        biotype,
        rank,
        hgvs_c,
        hgvs_p,
        cds_pos,
        aa_pos,
    ]

    key_material = "\x1f".join(
        [
            clean(chrom),
            clean(pos),
            clean(ref),
            clean(alt),
            effect,
            impact,
            gene,
            gene_id,
            product_name,
            strand,
            feature_id,
            biotype,
            rank,
            hgvs_c,
            hgvs_p,
            cds_pos,
            aa_pos,
        ]
    )
    site_key = hashlib.sha1(key_material.encode("utf-8")).hexdigest()
    return [site_key] + fields_without_key



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
