#!/usr/bin/env python3
"""Write TB-Profiler resistance results into a two-level Excel sheet."""

from pathlib import Path
import json
import argparse
import re
from typing import Dict, Any, List, Tuple, Optional

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
FIELDS = ["Pos", "GeneName", "Mutation", "Freq", "Confidence", "реф-й кодон", "кодон с заменой"]
DRUG_MAP = {d.lower(): d for d in DRUGS}
OTHER_VARIANTS_COL = ("", "Other Variants")  # MultiIndex tuple
DR_TYPE_COL = ("", "dr_type")
NT_SUB_RE = re.compile(r"^c\.(\d+)([ACGTN]+)>([ACGTN]+)$", re.IGNORECASE)
NT_DELINS_RE = re.compile(r"^c\.(\d+)_(\d+)del([ACGTN]+)ins([ACGTN]+)$", re.IGNORECASE)
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


class GeneModel:
    def __init__(self, key: str, strand: int, coding_seq: str, coding_positions: List[int]) -> None:
        self.key = key
        self.strand = strand
        self.coding_seq = coding_seq
        self.coding_positions = coding_positions

# ─────────── helpers ─────────────────────────────────────────────────────── #

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

        model = GeneModel(
            key=names[0],
            strand=int(feature.location.strand or 1),
            coding_seq=coding_seq,
            coding_positions=coding_positions,
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
                    codon_indices = [model.coding_positions.index(gp) for gp in codon_gpos]
                except ValueError:
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
        if model.coding_seq[start_idx : start_idx + len(ref)] != ref:
            return "", ""
    else:
        pos = var.get("pos")
        ref = str(var.get("ref") or "").upper()
        alt = str(var.get("alt") or "").upper()
        if not isinstance(pos, int) or not ref or not alt or len(ref) != len(alt):
            return "", ""

        variant_positions = list(range(pos, pos + len(ref)))
        try:
            variant_indices = [model.coding_positions.index(genome_pos) for genome_pos in variant_positions]
        except ValueError:
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

        if model.coding_seq[sorted_indices[0] : sorted_indices[0] + len(ref_seq)] != ref_seq:
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


def summarise(
    sample: Dict[str, Any],
    gene_index: Optional[Dict[str, GeneModel]] = None,
    block_name: str = "dr_variants",
    include_other_uncertain: bool = False,
) -> Dict[Tuple[str, str], str]:
    """Convert one TB‑Profiler JSON → dict suitable for a MultiIndex DataFrame.

    • Drug-specific data are stored under the selected top-level variant block.
    • "Other Variants" is simply the length of the top‑level
      ``other_variants`` list (0 if absent).
    • Frequency values are rounded to **two decimal places**.
    • Codon columns are calculated from the reference GenBank, following the
      same logic as tbprofiler_variant_table.py.
    """
    # Keep per-field lists so empty confidence values still occupy a slot.
    row_parts: Dict[Tuple[str, str], List[str]] = {
        (d, f): [] for d in DRUGS for f in FIELDS
    }

    # Pre‑populate the special column
    other_variants = str(len(sample.get("other_variants", [])))
    dr_type = str(sample.get("drtype") or "")

    # Helper to accumulate field values while optionally preserving blanks
    def acc(key: Tuple[str, str], val: Any, *, keep_blank: bool = False) -> None:
        if val is None:
            return
        text = "" if val == "" else str(val)
        if text == "" and not keep_blank:
            return
        row_parts[key].append(text)

    def combine_confidences(annotations: List[Dict[str, Any]]) -> str:
        """Match the TXT resistance-variants report: one variant entry per drug,
        with multiple labels for the same variant collapsed into one cell.
        Comments are only used when every annotation in the group has an empty
        confidence. If at least one confidence is present, comments from the
        blank-confidence annotations are ignored.
        """
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
        # Mirror the resistance-variants TXT block: one entry per variant/drug,
        # even if TB-Profiler attached multiple annotations for that drug.
            grouped_drugs = extract_annotations_for_variant(var, current_block)
            for drug_name, annotations in grouped_drugs.items():
                if current_block == "other_variants":
                    labels = [annotation_label(ann) for ann in annotations]
                    if not any("Uncertain significance" in label for label in labels):
                        continue
                d = DRUG_MAP.get(drug_name)
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
                acc((d, "Confidence"), combine_confidences(annotations), keep_blank=True)
                ref_codon, alt_codon = codons_from_variant(var, gene_index)
                acc((d, "реф-й кодон"), format_codons(ref_codon), keep_blank=True)
                acc((d, "кодон с заменой"), format_codons(alt_codon), keep_blank=True)

    row = {
        key: ("" if not any(v != "" for v in values) else ";".join(values))
        for key, values in row_parts.items()
    }
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


def write_excel(df: pd.DataFrame, out_path: Path) -> None:
    """Write the summary to Excel with a two‑level header.

    The first column is *Sample*, the second is *Other Variants*, followed by
    a block for every drug and its per-variant fields.
    """
    wb = xlsxwriter.Workbook(out_path)
    ws = wb.add_worksheet("TB‑Summary")

    hdr_fmt = wb.add_format({"bold": True, "align": "center", "valign": "vcenter", "border": 1})
    sub_fmt = wb.add_format({"align": "center", "valign": "vcenter", "border": 1})
    cell_fmt = wb.add_format({"border": 1})
    ref_codon_fmt = wb.add_format({"border": 1, "bg_color": "#DCE6D1"})
    alt_codon_fmt = wb.add_format({"border": 1, "bg_color": "#FFF2CC"})
    alt_codon_changed_fmt = wb.add_format({"border": 1, "bg_color": "#FFF2CC", "font_color": "#FF0000"})

    # ── header rows ──
    ws.write_blank(0, 0, "", hdr_fmt)        # (0,0) — top‑left corner
    ws.write(1, 0, "Sample", sub_fmt)         # row 1 header for Sample

    ws.write_blank(0, 1, "", hdr_fmt)        # blank top‑level for Other Variants
    ws.write(1, 1, "Other Variants", sub_fmt)
    ws.write_blank(0, 2, "", hdr_fmt)        # blank top-level for dr_type
    ws.write(1, 2, "dr_type", sub_fmt)

    col = 3  # start writing drug blocks from the fourth column
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
        ws.write(r, 2, row[DR_TYPE_COL], cell_fmt)
        col = 3
        for drug in DRUGS:
            for field in FIELDS:
                value = row[(drug, field)]
                if field == "реф-й кодон":
                    ws.write(r, col, value, ref_codon_fmt)
                elif field == "кодон с заменой":
                    ref_value = row[(drug, "реф-й кодон")]
                    segments = build_alt_codon_segments(value, ref_value, alt_codon_fmt, alt_codon_changed_fmt)
                    if segments:
                        ws.write_rich_string(r, col, *segments, alt_codon_fmt)
                    else:
                        ws.write(r, col, value, alt_codon_fmt)
                else:
                    ws.write(r, col, value, cell_fmt)
                col += 1

    # Keep Sample and Other Variants columns visible when scrolling
    ws.freeze_panes(2, 3)
    # Apply autofilter across the full width of the data
    ws.autofilter(1, 0, df.shape[0] + 1, df.shape[1])
    wb.close()


def write_excel_flat(df: pd.DataFrame, out_path: Path) -> None:
    wb = xlsxwriter.Workbook(out_path)
    ws = wb.add_worksheet("TB‑Summary")

    hdr_fmt = wb.add_format({"bold": True, "align": "center", "valign": "vcenter", "border": 1})
    cell_fmt = wb.add_format({"border": 1})
    ref_codon_fmt = wb.add_format({"border": 1, "bg_color": "#DCE6D1"})
    alt_codon_fmt = wb.add_format({"border": 1, "bg_color": "#FFF2CC"})
    alt_codon_changed_fmt = wb.add_format({"border": 1, "bg_color": "#FFF2CC", "font_color": "#FF0000"})

    ws.write(0, 0, "Sample", hdr_fmt)
    for col, column_name in enumerate(df.columns, start=1):
        ws.write(0, col, column_name, hdr_fmt)

    for r, (sample_id, row) in enumerate(df.iterrows(), start=1):
        ws.write(r, 0, sample_id, cell_fmt)
        for col, column_name in enumerate(df.columns, start=1):
            value = row[column_name]
            if column_name.endswith("_реф-й кодон"):
                ws.write(r, col, value, ref_codon_fmt)
            elif column_name.endswith("_кодон с заменой"):
                ref_name = column_name.replace("_кодон с заменой", "_реф-й кодон")
                ref_value = row.get(ref_name, "")
                segments = build_alt_codon_segments(value, ref_value, alt_codon_fmt, alt_codon_changed_fmt)
                if segments:
                    ws.write_rich_string(r, col, *segments, alt_codon_fmt)
                else:
                    ws.write(r, col, value, alt_codon_fmt)
            else:
                ws.write(r, col, value, cell_fmt)

    ws.freeze_panes(1, 3)
    ws.autofilter(0, 0, df.shape[0], df.shape[1])
    wb.close()


def flatten_columns(columns: pd.Index) -> List[str]:
    flat: List[str] = []
    for col in columns:
        if isinstance(col, tuple):
            if not col[0]:
                flat.append(col[1])
            else:
                flat.append(f"{col[0]}_{col[1]}")
        else:
            flat.append(str(col))
    return flat


def other_variants_output_path(out_path: Path) -> Path:
    return out_path.with_name(f"{out_path.stem}_other_variants{out_path.suffix}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="TB‑Profiler 2‑level Excel summary.")
    p.add_argument("-i", "--input", required=True, nargs="+", metavar="FILE.json", help="One or more TB‑Profiler JSON files.")
    p.add_argument("-o", "--output", metavar="FILE.xlsx", help="Output Excel file; leave blank to just print.")
    p.add_argument(
        "-g",
        "--genbank",
        metavar="FILE.gb",
        help="Reference GenBank file for codon calculation. Defaults to h37rv.gb in the current directory when present.",
    )
    p.add_argument(
        "--flat-header",
        action="store_true",
        help="Write a one-level table with columns like Isoniazid_GeneName instead of a two-level drug header.",
    )
    p.add_argument(
        "--also-other-variants",
        action="store_true",
        help="Additionally write the same table for the other_variants block with the suffix _other_variants.",
    )
    p.add_argument(
        "--include-other-uncertain",
        action="store_true",
        help="Also include Uncertain significance variants from other_variants in the main table.",
    )
    return p.parse_args()

# ─────────── main ────────────────────────────────────────────────────────── #

def main() -> None:
    args = parse_args()
    genbank_path = Path(args.genbank) if args.genbank else Path("h37rv.gb")
    gene_index = build_gene_index(genbank_path) if genbank_path.exists() else None

    def build_df(block_name: str) -> pd.DataFrame:
        rows: List[Dict[Tuple[str, str], str]] = []
        ids: List[str] = []
        for fname in args.input:
            p = Path(fname)
            sid = p.stem.split(".")[0]
            rows.append(
                summarise(
                    load_json(p),
                    gene_index=gene_index,
                    block_name=block_name,
                    include_other_uncertain=args.include_other_uncertain,
                )
            )
            ids.append(sid)

        columns = pd.MultiIndex.from_tuples([OTHER_VARIANTS_COL, DR_TYPE_COL]).append(
            pd.MultiIndex.from_product([DRUGS, FIELDS])
        )

        df = pd.DataFrame(rows, index=ids, columns=columns)
        df.index.name = "Sample"
        return df

    def emit_df(df: pd.DataFrame, out_path: Optional[Path]) -> None:
        if args.flat_header:
            flat = df.copy()
            flat.columns = flatten_columns(flat.columns)
            if out_path:
                write_excel_flat(flat, out_path)
                print(f"Saved {len(flat)} sample(s) ➜ {out_path}")
            else:
                print(flat.reset_index())
            return

        if out_path:
            write_excel(df, out_path)
            print(f"Saved {len(df)} sample(s) ➜ {out_path}")
        else:
            flat = df.copy()
            flat.columns = flatten_columns(flat.columns)
            print(flat.reset_index())

    main_df = build_df("dr_variants")
    main_out = Path(args.output) if args.output else None
    emit_df(main_df, main_out)

    if args.also_other_variants:
        other_df = build_df("other_variants")
        other_out = other_variants_output_path(Path(args.output)) if args.output else None
        emit_df(other_df, other_out)


if __name__ == "__main__":
    main()
