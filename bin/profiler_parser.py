#!/usr/bin/env python3
"""Write TB-Profiler resistance results into Excel."""

from __future__ import annotations

from pathlib import Path
from glob import glob
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
import argparse
import json
import os
import re
import sys
import time
from typing import Dict, Any, List, Tuple, Optional, Iterable

import pandas as pd
import xlsxwriter

# ─────────── constant data ──────────────────────────────────────────────── #
DRUGS = [
    "Rifampicin", "Isoniazid", "Ethambutol", "Pyrazinamide",
    "Moxifloxacin", "Levofloxacin", "Bedaquiline", "Delamanid",
    "Pretomanid", "Linezolid", "Streptomycin", "Amikacin",
    "Kanamycin", "Capreomycin", "Clofazimine", "Ethionamide",
    "Para-aminosalicylic_acid", "Cycloserine",
]
DRUG_MAP = {d.lower(): d for d in DRUGS}
OTHER_VARIANTS_COL = ("", "Other Variants")
DR_TYPE_COL = ("", "dr_type")
RESISTANCE_FLAG_COL = "Resistance"
RESISTANCE_CONFIDENCE_LABELS = {"Assoc w R", "Assoc w R - Interim"}
NT_SUB_RE = re.compile(r"^c\.(\d+)([ACGTN]+)>([ACGTN]+)$", re.IGNORECASE)
NT_DELINS_RE = re.compile(r"^c\.(\d+)_(\d+)del([ACGTN]+)ins([ACGTN]+)$", re.IGNORECASE)
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

FIELDS = [RESISTANCE_FLAG_COL, "Pos", "GeneName", "Mutation", "Freq", "Confidence", "реф-й кодон", "кодон с заменой"]

ORDERED_COLUMNS: List[Tuple[str, str]] = [OTHER_VARIANTS_COL, DR_TYPE_COL] + [
    (drug, field) for drug in DRUGS for field in FIELDS
]


# ─────────── worker globals ─────────────────────────────────────────────── #
_WORKER_GENE_INDEX: Optional[Dict[str, "GeneModel"]] = None
_WORKER_INCLUDE_OTHER_UNCERTAIN = False
_WORKER_RETURN_OTHER = False


@dataclass
class GeneModel:
    key: str
    strand: int
    coding_seq: str
    coding_positions: List[int]
    pos_to_index: Dict[int, int]


# ─────────── helpers ────────────────────────────────────────────────────── #

def load_json(path: Path) -> Dict[str, Any]:
    with path.open() as fh:
        return json.load(fh)


def feature_names(feature: Any) -> List[str]:
    names: List[str] = []
    qualifiers = feature.qualifiers
    for key in ("gene", "locus_tag", "old_locus_tag", "gene_synonym", "label"):
        for value in qualifiers.get(key, []):
            if value and value not in names:
                names.append(value)
    return names


def coding_positions_for_feature(feature: Any) -> List[int]:
    loc = feature.location
    parts = list(getattr(loc, "parts", [loc]))
    strand = int(loc.strand or 1)
    ordered_parts = parts if strand == 1 else list(reversed(parts))

    positions: List[int] = []
    for part in ordered_parts:
        start = int(part.start)
        end = int(part.end)
        if strand == 1:
            positions.extend(range(start + 1, end + 1))
        else:
            positions.extend(range(end, start, -1))
    return positions


def build_gene_index(genbank_path: Path) -> Dict[str, GeneModel]:
    try:
        from Bio import SeqIO
    except ImportError as exc:
        raise SystemExit(
            "Biopython is required for codon calculation from GenBank. "
            "Install the 'biopython' package or run without a GenBank file."
        ) from exc

    record = SeqIO.read(str(genbank_path), "genbank")
    gene_index: Dict[str, GeneModel] = {}

    for feature in record.features:
        if feature.type not in {"CDS", "gene"}:
            continue

        names = feature_names(feature)
        if not names:
            continue

        coding_seq = str(feature.extract(record.seq)).upper()
        coding_positions = coding_positions_for_feature(feature)
        if not coding_seq or len(coding_seq) != len(coding_positions):
            continue

        pos_to_index = {pos: i for i, pos in enumerate(coding_positions)}
        model = GeneModel(
            key=names[0],
            strand=int(feature.location.strand or 1),
            coding_seq=coding_seq,
            coding_positions=coding_positions,
            pos_to_index=pos_to_index,
        )
        for name in names:
            gene_index[name] = model

    return gene_index


def complement_base(base: str) -> str:
    return base.translate(COMPLEMENT)


def gene_model_for_variant(var: Dict[str, Any], gene_index: Dict[str, GeneModel]) -> Optional[GeneModel]:
    for key in (var.get("gene_id"), var.get("gene_name"), var.get("feature_id"), var.get("locus_tag")):
        if key and key in gene_index:
            return gene_index[key]
    return None


def coding_change_from_nucleotide_change(var: Dict[str, Any]) -> Optional[Tuple[int, str, str]]:
    nucleotide_change = str(var.get("nucleotide_change") or "").strip()
    if not nucleotide_change:
        return None

    m = NT_SUB_RE.match(nucleotide_change)
    if m:
        start, ref, alt = m.groups()
        return int(start) - 1, ref.upper(), alt.upper()

    m = NT_DELINS_RE.match(nucleotide_change)
    if m:
        start, end, ref, alt = m.groups()
        if (int(end) - int(start) + 1) != len(ref):
            return None
        return int(start) - 1, ref.upper(), alt.upper()

    return None


def codons_from_variant(var: Dict[str, Any], gene_index: Optional[Dict[str, GeneModel]]) -> Tuple[str, str]:
    if gene_index is None:
        return "", ""

    model = gene_model_for_variant(var, gene_index)
    if model is None:
        return "", ""

    coding_change = coding_change_from_nucleotide_change(var)
    if coding_change is not None:
        start_idx, ref, alt = coding_change
        if len(ref) != len(alt) or start_idx < 0:
            return "", ""

        sorted_indices = list(range(start_idx, start_idx + len(ref)))

        genome_var_pos = var.get("pos")
        if isinstance(genome_var_pos, int):
            anchor_idx = sorted_indices[-1] if model.strand == -1 else sorted_indices[0]
            if anchor_idx < len(model.coding_positions) and model.coding_positions[anchor_idx] != genome_var_pos:
                codon_start_tbp = (sorted_indices[0] // 3) * 3
                codon_end_tbp = ((sorted_indices[-1] // 3) + 1) * 3

                if model.strand == -1:
                    tbp_gene_end = genome_var_pos + sorted_indices[-1]
                    codon_gpos = [tbp_gene_end - i for i in range(codon_start_tbp, codon_end_tbp)]
                else:
                    tbp_gene_start = genome_var_pos - sorted_indices[0]
                    codon_gpos = [tbp_gene_start + i for i in range(codon_start_tbp, codon_end_tbp)]

                try:
                    codon_indices = [model.pos_to_index[gp] for gp in codon_gpos]
                except KeyError:
                    return "", ""

                ref_codon = "".join(model.coding_seq[i] for i in codon_indices)
                alt_by_tbp = dict(zip(sorted_indices, alt))
                alt_codon_list = list(ref_codon)
                for j, tbp_idx in enumerate(range(codon_start_tbp, codon_end_tbp)):
                    if tbp_idx in alt_by_tbp:
                        alt_codon_list[j] = alt_by_tbp[tbp_idx]
                return ref_codon, "".join(alt_codon_list)

        if sorted_indices[-1] >= len(model.coding_seq):
            return "", ""
        if model.coding_seq[start_idx:start_idx + len(ref)] != ref:
            return "", ""
    else:
        pos = var.get("pos")
        ref = str(var.get("ref") or "").upper()
        alt = str(var.get("alt") or "").upper()
        if not isinstance(pos, int) or not ref or not alt or len(ref) != len(alt):
            return "", ""

        variant_positions = list(range(pos, pos + len(ref)))
        try:
            variant_indices = [model.pos_to_index[genome_pos] for genome_pos in variant_positions]
        except KeyError:
            return "", ""

        sorted_indices = sorted(variant_indices)
        if sorted_indices != list(range(sorted_indices[0], sorted_indices[0] + len(ref))):
            return "", ""

        if model.strand == 1:
            ref_seq = ref
            alt_seq = alt
        else:
            ref_seq = complement_base(ref[::-1])
            alt_seq = complement_base(alt[::-1])

        if model.coding_seq[sorted_indices[0]:sorted_indices[0] + len(ref_seq)] != ref_seq:
            return "", ""

        ref = ref_seq
        alt = alt_seq

    codon_start = (sorted_indices[0] // 3) * 3
    codon_end = ((sorted_indices[-1] // 3) + 1) * 3
    if codon_end > len(model.coding_seq):
        return "", ""

    ref_codon = model.coding_seq[codon_start:codon_end]
    alt_by_index = dict(zip(sorted_indices, alt))
    alt_codon_list = list(ref_codon)

    for idx in sorted_indices:
        alt_codon_list[idx - codon_start] = alt_by_index[idx]

    return ref_codon, "".join(alt_codon_list)


def format_codons(codon: str) -> str:
    return " ".join(codon[i:i + 3] for i in range(0, len(codon), 3))


def extract_annotations_for_variant(var: Dict[str, Any], block_name: str) -> Dict[str, List[Dict[str, Any]]]:
    grouped: Dict[str, List[Dict[str, Any]]] = {}

    if block_name == "dr_variants":
        raw_items = var.get("drugs", [])
    else:
        raw_items = var.get("annotation", [])

    if isinstance(raw_items, list):
        for item in raw_items:
            if not isinstance(item, dict):
                continue
            drug_name = str(item.get("drug") or "").strip().lower()
            if not drug_name:
                continue
            grouped.setdefault(drug_name, []).append(item)

    if not grouped:
        for drug_name in var.get("gene_associated_drugs") or []:
            drug_key = str(drug_name).strip().lower()
            if not drug_key:
                continue
            grouped.setdefault(drug_key, []).append({"drug": drug_key, "confidence": ""})

    return grouped


def annotation_label(annotation: Dict[str, Any]) -> str:
    confidence = str(annotation.get("confidence") or "").strip()
    if confidence:
        return confidence
    return str(annotation.get("comment") or "").strip()


def is_resistance_annotation(annotation: Dict[str, Any]) -> bool:
    return annotation_label(annotation) in RESISTANCE_CONFIDENCE_LABELS


def summarise(
    sample: Dict[str, Any],
    gene_index: Optional[Dict[str, GeneModel]] = None,
    block_name: str = "dr_variants",
    include_other_uncertain: bool = False,
) -> Dict[Tuple[str, str], str]:
    row_parts: Dict[Tuple[str, str], List[str]] = {
        (d, f): [] for d in DRUGS for f in FIELDS
    }
    resistance_flags: Dict[str, bool] = {drug: False for drug in DRUGS}

    other_variants = str(len(sample.get("other_variants", [])))
    dr_type = str(sample.get("drtype") or "")

    def acc(key: Tuple[str, str], val: Any, *, keep_blank: bool = False) -> None:
        if val is None:
            return
        text = "" if val == "" else str(val)
        if text == "" and not keep_blank:
            return
        row_parts[key].append(text)

    def combine_confidences(annotations: List[Dict[str, Any]]) -> str:
        seen = set()
        confidence_values: List[str] = []
        comment_values: List[str] = []

        for ann in annotations:
            confidence = str(ann.get("confidence") or "").strip()
            comment = str(ann.get("comment") or "").strip()

            if confidence:
                if confidence not in seen:
                    seen.add(confidence)
                    confidence_values.append(confidence)
                continue

            if comment and comment not in seen:
                seen.add(comment)
                comment_values.append(comment)

        if confidence_values:
            return "|".join(confidence_values)
        if comment_values:
            return "|".join(comment_values)
        return ""

    variant_blocks = [block_name]
    if block_name == "dr_variants" and include_other_uncertain:
        variant_blocks.append("other_variants")

    for current_block in variant_blocks:
        for var in sample.get(current_block, []):
            grouped_drugs = extract_annotations_for_variant(var, current_block)
            for drug_name, annotations in grouped_drugs.items():
                if current_block == "other_variants":
                    labels = [annotation_label(ann) for ann in annotations]
                    if not any("Uncertain significance" in label for label in labels):
                        continue

                d = DRUG_MAP.get(drug_name)
                if not d:
                    continue

                if any(is_resistance_annotation(annotation) for annotation in annotations):
                    resistance_flags[d] = True

                freq_val = var.get("freq")
                if freq_val not in (None, ""):
                    try:
                        freq_val = f"{float(freq_val):.2f}"
                    except (ValueError, TypeError):
                        pass

                acc((d, "Pos"), var.get("pos"))
                acc((d, "GeneName"), var.get("gene_name"))
                acc((d, "Mutation"), var.get("change"))
                acc((d, "Freq"), freq_val)
                acc((d, "Confidence"), combine_confidences(annotations), keep_blank=True)

                ref_codon, alt_codon = codons_from_variant(var, gene_index)
                acc((d, "реф-й кодон"), format_codons(ref_codon), keep_blank=True)
                acc((d, "кодон с заменой"), format_codons(alt_codon), keep_blank=True)

    row = {
        key: ("" if not any(v != "" for v in values) else ";".join(values))
        for key, values in row_parts.items()
    }
    for drug, has_resistance in resistance_flags.items():
        row[(drug, RESISTANCE_FLAG_COL)] = "R" if has_resistance else ""
    row[OTHER_VARIANTS_COL] = other_variants
    row[DR_TYPE_COL] = dr_type
    return row


def build_alt_codon_segments(
    alt_value: str,
    ref_value: str,
    default_fmt: Any,
    changed_fmt: Any,
) -> Optional[List[Any]]:
    if not alt_value or not ref_value or len(alt_value) != len(ref_value):
        return None

    segments: List[Any] = []
    current_fmt = changed_fmt if alt_value[0] != ref_value[0] else default_fmt
    current_text = alt_value[0]

    for i in range(1, len(alt_value)):
        fmt = changed_fmt if alt_value[i] != ref_value[i] else default_fmt
        if fmt is current_fmt:
            current_text += alt_value[i]
        else:
            segments.extend([current_fmt, current_text])
            current_fmt = fmt
            current_text = alt_value[i]

    segments.extend([current_fmt, current_text])
    if len(segments) < 4:
        return None
    return segments


def flatten_columns(columns: Iterable[Tuple[str, str]]) -> List[str]:
    flat: List[str] = []
    for col in columns:
        if not col[0]:
            flat.append(col[1])
        else:
            flat.append(f"{col[0]}_{col[1]}")
    return flat


FLAT_COLUMN_NAMES = flatten_columns(ORDERED_COLUMNS)


def other_variants_output_path(out_path: Path) -> Path:
    return out_path.with_name(f"{out_path.stem}_other_variants{out_path.suffix}")


def read_input_list(path: Path) -> List[str]:
    items: List[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            items.append(line)
    return items


def resolve_inputs(input_items: List[str], input_list: Optional[Path] = None) -> List[Path]:
    raw_items = list(input_items or [])
    if input_list is not None:
        raw_items.extend(read_input_list(input_list))

    resolved: List[Path] = []
    seen = set()

    for item in raw_items:
        item = item.strip()
        if not item:
            continue

        matches: List[Path]
        p = Path(item)

        if p.is_dir():
            matches = sorted(p.rglob("*.results.json"))
        elif any(ch in item for ch in "*?[]"):
            matches = [Path(x) for x in glob(item, recursive=True)]
        else:
            matches = [p]

        for match in matches:
            if match.is_file():
                key = str(match.resolve())
                if key not in seen:
                    seen.add(key)
                    resolved.append(match)

    if not resolved:
        raise SystemExit("No input JSON files found.")

    return sorted(resolved)


def init_worker(
    gene_index: Optional[Dict[str, GeneModel]],
    include_other_uncertain: bool,
    return_other: bool,
) -> None:
    global _WORKER_GENE_INDEX, _WORKER_INCLUDE_OTHER_UNCERTAIN, _WORKER_RETURN_OTHER
    _WORKER_GENE_INDEX = gene_index
    _WORKER_INCLUDE_OTHER_UNCERTAIN = include_other_uncertain
    _WORKER_RETURN_OTHER = return_other


def process_one_json(path_str: str) -> Tuple[str, Dict[Tuple[str, str], str], Optional[Dict[Tuple[str, str], str]]]:
    path = Path(path_str)
    try:
        sample = load_json(path)
        sample_id = path.stem.split(".")[0]

        main_row = summarise(
            sample,
            gene_index=_WORKER_GENE_INDEX,
            block_name="dr_variants",
            include_other_uncertain=_WORKER_INCLUDE_OTHER_UNCERTAIN,
        )

        other_row = None
        if _WORKER_RETURN_OTHER:
            other_row = summarise(
                sample,
                gene_index=_WORKER_GENE_INDEX,
                block_name="other_variants",
                include_other_uncertain=False,
            )

        return sample_id, main_row, other_row
    except Exception as exc:
        raise RuntimeError(f"Failed to process {path}: {exc}") from exc


def iter_processed(
    input_files: List[Path],
    gene_index: Optional[Dict[str, GeneModel]],
    include_other_uncertain: bool,
    return_other: bool,
    jobs: int,
    chunksize: int,
) -> Iterable[Tuple[str, Dict[Tuple[str, str], str], Optional[Dict[Tuple[str, str], str]]]]:
    input_strings = [str(p) for p in input_files]

    if jobs <= 1:
        init_worker(gene_index, include_other_uncertain, return_other)
        for item in input_strings:
            yield process_one_json(item)
        return

    with ProcessPoolExecutor(
        max_workers=jobs,
        initializer=init_worker,
        initargs=(gene_index, include_other_uncertain, return_other),
    ) as executor:
        yield from executor.map(process_one_json, input_strings, chunksize=chunksize)


def iter_with_progress(
    iterator: Iterable[Any],
    total: int,
    enabled: bool = True,
    label: str = "Processing",
) -> Iterable[Any]:
    if not enabled:
        yield from iterator
        return

    start = time.monotonic()
    last_print = start

    for idx, item in enumerate(iterator, start=1):
        now = time.monotonic()
        if idx == 1 or idx == total or (now - last_print) >= 1.0:
            rate = idx / max(now - start, 1e-9)
            pct = idx * 100.0 / total
            print(
                f"\r{label}: {idx}/{total} ({pct:5.1f}%) | {rate:,.1f} files/s",
                end="",
                file=sys.stderr,
                flush=True,
            )
            last_print = now
        yield item

    print(file=sys.stderr, flush=True)


class SummaryWorkbookWriter:
    def __init__(self, out_path: Path, flat_header: bool) -> None:
        self.out_path = out_path
        self.flat_header = flat_header
        self.workbook = xlsxwriter.Workbook(str(out_path))
        self.ws = self.workbook.add_worksheet("TB-Summary")

        self.hdr_fmt = self.workbook.add_format({"bold": True, "align": "center", "valign": "vcenter", "border": 1})
        self.sub_fmt = self.workbook.add_format({"align": "center", "valign": "vcenter", "border": 1})
        self.cell_fmt = self.workbook.add_format({"border": 1})
        self.ref_codon_fmt = self.workbook.add_format({"border": 1, "bg_color": "#DCE6D1"})
        self.alt_codon_fmt = self.workbook.add_format({"border": 1, "bg_color": "#FFF2CC"})
        self.alt_codon_changed_fmt = self.workbook.add_format({"border": 1, "bg_color": "#FFF2CC", "font_color": "#FF0000"})

        if self.flat_header:
            self.data_start_row = 1
            self._write_flat_header()
        else:
            self.data_start_row = 2
            self._write_multi_header()

    def _write_multi_header(self) -> None:
        self.ws.write_blank(0, 0, "", self.hdr_fmt)
        self.ws.write(1, 0, "Sample", self.sub_fmt)

        self.ws.write_blank(0, 1, "", self.hdr_fmt)
        self.ws.write(1, 1, "Other Variants", self.sub_fmt)

        self.ws.write_blank(0, 2, "", self.hdr_fmt)
        self.ws.write(1, 2, "dr_type", self.sub_fmt)

        col = 3
        for drug in DRUGS:
            start, end = col, col + len(FIELDS) - 1
            self.ws.merge_range(0, start, 0, end, drug, self.hdr_fmt)
            for i, field in enumerate(FIELDS):
                self.ws.write(1, start + i, field, self.sub_fmt)
            col = end + 1

    def _write_flat_header(self) -> None:
        self.ws.write(0, 0, "Sample", self.hdr_fmt)
        for col, column_name in enumerate(FLAT_COLUMN_NAMES, start=1):
            self.ws.write(0, col, column_name, self.hdr_fmt)

    def write_record(self, row_number: int, sample_id: str, row: Dict[Tuple[str, str], str]) -> None:
        excel_row = self.data_start_row + row_number - 1

        self.ws.write(excel_row, 0, sample_id, self.cell_fmt)

        if self.flat_header:
            for col, key in enumerate(ORDERED_COLUMNS, start=1):
                value = row[key]
                if key[1] == "реф-й кодон":
                    self.ws.write(excel_row, col, value, self.ref_codon_fmt)
                elif key[1] == "кодон с заменой":
                    ref_key = (key[0], "реф-й кодон")
                    ref_value = row.get(ref_key, "")
                    segments = build_alt_codon_segments(value, ref_value, self.alt_codon_fmt, self.alt_codon_changed_fmt)
                    if segments:
                        self.ws.write_rich_string(excel_row, col, *segments, self.alt_codon_fmt)
                    else:
                        self.ws.write(excel_row, col, value, self.alt_codon_fmt)
                else:
                    self.ws.write(excel_row, col, value, self.cell_fmt)
            return

        self.ws.write(excel_row, 1, row[OTHER_VARIANTS_COL], self.cell_fmt)
        self.ws.write(excel_row, 2, row[DR_TYPE_COL], self.cell_fmt)

        col = 3
        for drug in DRUGS:
            for field in FIELDS:
                key = (drug, field)
                value = row[key]
                if field == "реф-й кодон":
                    self.ws.write(excel_row, col, value, self.ref_codon_fmt)
                elif field == "кодон с заменой":
                    ref_value = row[(drug, "реф-й кодон")]
                    segments = build_alt_codon_segments(value, ref_value, self.alt_codon_fmt, self.alt_codon_changed_fmt)
                    if segments:
                        self.ws.write_rich_string(excel_row, col, *segments, self.alt_codon_fmt)
                    else:
                        self.ws.write(excel_row, col, value, self.alt_codon_fmt)
                else:
                    self.ws.write(excel_row, col, value, self.cell_fmt)
                col += 1

    def close(self, total_rows: int) -> None:
        if self.flat_header:
            self.ws.freeze_panes(1, 3)
            self.ws.autofilter(0, 0, total_rows, len(ORDERED_COLUMNS))
        else:
            self.ws.freeze_panes(2, 3)
            self.ws.autofilter(1, 0, total_rows + 1, len(ORDERED_COLUMNS))
        self.workbook.close()


def records_to_dataframe(
    records: List[Tuple[str, Dict[Tuple[str, str], str]]]
) -> pd.DataFrame:
    ids = [sample_id for sample_id, _row in records]
    rows = [row for _sample_id, row in records]
    df = pd.DataFrame(rows, index=ids, columns=pd.MultiIndex.from_tuples(ORDERED_COLUMNS))
    df.index.name = "Sample"
    return df


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="TB-Profiler Excel summary.")
    p.add_argument(
        "-i", "--input",
        nargs="*",
        default=[],
        metavar="INPUT",
        help="Input file(s), directory(ies), or quoted glob(s). "
             "If a directory is given, *.results.json are searched recursively.",
    )
    p.add_argument(
        "--input-list",
        metavar="FILE.txt",
        help="Text file with one input file/directory/glob per line.",
    )
    p.add_argument("-o", "--output", metavar="FILE.xlsx", help="Output Excel file; leave blank to print.")
    p.add_argument(
        "-g", "--genbank",
        metavar="FILE.gb",
        help="Reference GenBank file for codon calculation. "
             "Defaults to h37rv.gb in current directory when present.",
    )
    p.add_argument(
        "--flat-header",
        action="store_true",
        help="Write one-level header like Isoniazid_GeneName.",
    )
    p.add_argument(
        "--also-other-variants",
        action="store_true",
        help="Additionally write the other_variants table with suffix _other_variants.",
    )
    p.add_argument(
        "--include-other-uncertain",
        action="store_true",
        help="Also include Uncertain significance variants from other_variants in main table.",
    )
    p.add_argument(
        "-j", "--jobs",
        type=int,
        default=max(1, min(8, os.cpu_count() or 1)),
        help="Number of worker processes. Default: min(8, CPU count).",
    )
    p.add_argument(
        "--chunksize",
        type=int,
        default=100,
        help="Chunk size for multiprocessing. Default: 100.",
    )
    p.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress output to stderr.",
    )
    args = p.parse_args()

    if not args.input and not args.input_list:
        p.error("provide --input and/or --input-list")

    if args.jobs < 1:
        p.error("--jobs must be >= 1")
    if args.chunksize < 1:
        p.error("--chunksize must be >= 1")

    return args


def main() -> None:
    args = parse_args()

    genbank_path = Path(args.genbank) if args.genbank else Path("h37rv.gb")
    gene_index: Optional[Dict[str, GeneModel]] = None

    if args.genbank:
        if not genbank_path.is_file():
            raise SystemExit(f"GenBank file not found: {genbank_path}")
        gene_index = build_gene_index(genbank_path)
    else:
        if genbank_path.is_file():
            gene_index = build_gene_index(genbank_path)

    input_list_path = Path(args.input_list) if args.input_list else None
    input_files = resolve_inputs(args.input, input_list_path)
    total = len(input_files)

    processed_iter = iter_processed(
        input_files=input_files,
        gene_index=gene_index,
        include_other_uncertain=args.include_other_uncertain,
        return_other=args.also_other_variants,
        jobs=args.jobs,
        chunksize=args.chunksize,
    )
    processed_iter = iter_with_progress(
        processed_iter,
        total=total,
        enabled=not args.no_progress,
        label="Processing",
    )

    if args.output:
        main_out = Path(args.output)
        other_out = other_variants_output_path(main_out) if args.also_other_variants else None

        main_writer = SummaryWorkbookWriter(main_out, flat_header=args.flat_header)
        other_writer = SummaryWorkbookWriter(other_out, flat_header=args.flat_header) if other_out else None

        count = 0
        try:
            for count, (sample_id, main_row, other_row) in enumerate(processed_iter, start=1):
                main_writer.write_record(count, sample_id, main_row)
                if other_writer is not None and other_row is not None:
                    other_writer.write_record(count, sample_id, other_row)
        finally:
            main_writer.close(count)
            if other_writer is not None:
                other_writer.close(count)

        print(f"Saved {count} sample(s) ➜ {main_out}")
        if other_out is not None:
            print(f"Saved {count} sample(s) ➜ {other_out}")
        return

    main_records: List[Tuple[str, Dict[Tuple[str, str], str]]] = []
    other_records: List[Tuple[str, Dict[Tuple[str, str], str]]] = []

    for sample_id, main_row, other_row in processed_iter:
        main_records.append((sample_id, main_row))
        if other_row is not None:
            other_records.append((sample_id, other_row))

    main_df = records_to_dataframe(main_records)
    flat_main = main_df.copy()
    flat_main.columns = flatten_columns(ORDERED_COLUMNS)
    print(flat_main.reset_index())

    if args.also_other_variants:
        other_df = records_to_dataframe(other_records)
        flat_other = other_df.copy()
        flat_other.columns = flatten_columns(ORDERED_COLUMNS)
        print("\n=== other_variants ===")
        print(flat_other.reset_index())


if __name__ == "__main__":
    main()
