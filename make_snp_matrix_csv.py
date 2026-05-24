#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
from pathlib import Path

VCF_EXTS = ('.vcf', '.vcf.gz')
VALID_SAMPLE_RE = re.compile(r'^[A-Za-z0-9._-]+$')


def is_vcf(path: Path) -> bool:
    name = path.name.lower()
    return any(name.endswith(ext) for ext in VCF_EXTS)


def sample_id_from_vcf(path: Path) -> str:
    name = path.name
    if name.lower().endswith('.vcf.gz'):
        sample = name[:-7]
    elif name.lower().endswith('.vcf'):
        sample = name[:-4]
    else:
        sample = path.stem

    parts = sample.split('__')
    if len(parts) == 2 and parts[0] == parts[1]:
        return parts[0]
    return sample


def read_vcf_sample_id(path: Path) -> str:
    opener = gzip.open if path.name.lower().endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('#CHROM\t'):
                fields = line.rstrip('\n').split('\t')
                samples = fields[9:]
                if len(samples) != 1:
                    raise SystemExit(
                        f"ERROR: '{path}' must contain exactly one sample column, found {len(samples)}"
                    )
                return samples[0]
    raise SystemExit(f"ERROR: VCF header line '#CHROM' not found in '{path}'")


def sanitize_sample_id(sample: str) -> str:
    sample = re.sub(r'[^A-Za-z0-9._-]+', '_', sample.strip())
    sample = sample.strip('._-')
    return sample


def find_vcfs(in_dir: Path, recursive: bool) -> list[Path]:
    iterator = in_dir.rglob('*') if recursive else in_dir.glob('*')
    return sorted(path.resolve() for path in iterator if path.is_file() and is_vcf(path))


def build_rows(in_dirs: list[Path], recursive: bool, sanitize: bool, sample_source: str) -> list[list[str]]:
    files = []
    for in_dir in in_dirs:
        dir_files = find_vcfs(in_dir, recursive)
        if not dir_files:
            raise SystemExit(f"ERROR: VCF files (*.vcf or *.vcf.gz) not found in '{in_dir}'")
        files.extend(dir_files)

    if not files:
        raise SystemExit("ERROR: VCF files (*.vcf or *.vcf.gz) not found in input directories")

    rows = []
    seen_samples = {}

    for vcf in files:
        sample = read_vcf_sample_id(vcf) if sample_source == 'header' else sample_id_from_vcf(vcf)
        if sanitize:
            sample = sanitize_sample_id(sample)

        if not sample:
            raise SystemExit(f"ERROR: cannot derive sample ID from '{vcf}'")
        if not VALID_SAMPLE_RE.fullmatch(sample):
            raise SystemExit(
                f"ERROR: sample ID '{sample}' from '{vcf.name}' contains unsupported characters. "
                "Rename the file or rerun with --sanitize-sample-ids."
            )
        if sample in seen_samples:
            raise SystemExit(
                f"ERROR: duplicate sample ID '{sample}' derived from both "
                f"'{seen_samples[sample]}' and '{vcf}'"
            )

        seen_samples[sample] = vcf
        rows.append([sample, str(vcf)])

    rows.sort(key=lambda row: row[0])
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(
        description='Create sample,vcf CSV input for snp_matrix.nf from one or more directories of VCF files.'
    )
    ap.add_argument(
        '-i',
        '--input',
        nargs='+',
        required=True,
        help='Input directory or directories with *.vcf or *.vcf.gz files',
    )
    ap.add_argument('-o', '--output', required=True, help='Output CSV path')
    ap.add_argument(
        '--no-recursive',
        action='store_true',
        help='Only scan the input directory itself; by default subdirectories are scanned too',
    )
    ap.add_argument(
        '--sanitize-sample-ids',
        action='store_true',
        help="Replace characters unsupported by snp_matrix.nf sample IDs with '_'",
    )
    ap.add_argument(
        '--sample-source',
        choices=('header', 'filename'),
        default='header',
        help='Where to get sample IDs from. Default: VCF header sample column',
    )
    args = ap.parse_args()

    in_dirs = [Path(path).resolve() for path in args.input]
    for in_dir in in_dirs:
        if not in_dir.is_dir():
            raise SystemExit(f"ERROR: '{in_dir}' is not a directory")

    rows = build_rows(
        in_dirs=in_dirs,
        recursive=not args.no_recursive,
        sanitize=args.sanitize_sample_ids,
        sample_source=args.sample_source,
    )

    out_path = Path(args.output).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open('w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['sample', 'vcf'])
        writer.writerows(rows)

    print(f"Done: {out_path} (samples: {len(rows)})")


if __name__ == '__main__':
    main()
