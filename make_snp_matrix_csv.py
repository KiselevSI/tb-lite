#!/usr/bin/env python3
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import gzip
import re
import sys
from pathlib import Path

VCF_EXTS = ('.vcf', '.vcf.gz')
VALID_SAMPLE_RE = re.compile(r'^[A-Za-z0-9._-]+$')
PROGRESS_INTERVAL = 1000


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
                    raise ValueError(
                        f"ERROR: '{path}' must contain exactly one sample column, found {len(samples)}"
                    )
                return samples[0]
    raise ValueError(f"ERROR: VCF header line '#CHROM' not found in '{path}'")


def sanitize_sample_id(sample: str) -> str:
    sample = re.sub(r'[^A-Za-z0-9._-]+', '_', sample.strip())
    sample = sample.strip('._-')
    return sample


def find_vcfs(in_dir: Path, recursive: bool) -> list[Path]:
    iterator = in_dir.rglob('*') if recursive else in_dir.glob('*')
    files = []
    for path in iterator:
        if path.is_file() and is_vcf(path):
            files.append(path.resolve())
            if len(files) % PROGRESS_INTERVAL == 0:
                print(f"Scanning {in_dir}: found {len(files)} VCF files", file=sys.stderr)
    print(f"Scanning {in_dir}: found {len(files)} VCF files", file=sys.stderr)
    return sorted(files)


def report_progress(done: int, total: int) -> None:
    if done == total or done % PROGRESS_INTERVAL == 0:
        print(f"Processed {done}/{total} VCF files", file=sys.stderr)


def build_row(vcf: Path, sanitize: bool, sample_source: str) -> list[str]:
    sample = read_vcf_sample_id(vcf) if sample_source == 'header' else sample_id_from_vcf(vcf)
    if sanitize:
        sample = sanitize_sample_id(sample)

    if not sample:
        raise ValueError(f"ERROR: cannot derive sample ID from '{vcf}'")
    if not VALID_SAMPLE_RE.fullmatch(sample):
        raise ValueError(
            f"ERROR: sample ID '{sample}' from '{vcf.name}' contains unsupported characters. "
            "Rename the file or rerun with --sanitize-sample-ids."
        )
    return [sample, str(vcf)]


def validate_unique_samples(rows: list[list[str]]) -> list[list[str]]:
    rows.sort(key=lambda row: (row[0], row[1]))
    seen_samples = {}

    for sample, vcf in rows:
        if sample in seen_samples:
            raise SystemExit(
                f"ERROR: duplicate sample ID '{sample}' derived from both "
                f"'{seen_samples[sample]}' and '{vcf}'"
            )
        seen_samples[sample] = vcf

    return rows


def build_rows(in_dirs: list[Path], recursive: bool, sanitize: bool, sample_source: str, jobs: int) -> list[list[str]]:
    files = []
    for in_dir in in_dirs:
        dir_files = find_vcfs(in_dir, recursive)
        if not dir_files:
            raise SystemExit(f"ERROR: VCF files (*.vcf or *.vcf.gz) not found in '{in_dir}'")
        files.extend(dir_files)

    if not files:
        raise SystemExit("ERROR: VCF files (*.vcf or *.vcf.gz) not found in input directories")

    print(f"Found {len(files)} VCF files", file=sys.stderr)
    rows = []

    if sample_source == 'header' and jobs > 1:
        with ThreadPoolExecutor(max_workers=jobs) as executor:
            futures = [executor.submit(build_row, vcf, sanitize, sample_source) for vcf in files]
            for done, future in enumerate(as_completed(futures), start=1):
                try:
                    rows.append(future.result())
                except ValueError as exc:
                    raise SystemExit(str(exc)) from exc
                report_progress(done, len(files))
    else:
        for done, vcf in enumerate(files, start=1):
            try:
                rows.append(build_row(vcf, sanitize, sample_source))
            except ValueError as exc:
                raise SystemExit(str(exc)) from exc
            report_progress(done, len(files))

    return validate_unique_samples(rows)


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
    ap.add_argument(
        '--jobs',
        type=int,
        default=1,
        help='Number of parallel workers for --sample-source header. Ignored by --sample-source filename. Default: 1',
    )
    args = ap.parse_args()

    if args.jobs < 1:
        raise SystemExit("ERROR: --jobs must be a positive integer")

    in_dirs = [Path(path).resolve() for path in args.input]
    for in_dir in in_dirs:
        if not in_dir.is_dir():
            raise SystemExit(f"ERROR: '{in_dir}' is not a directory")

    rows = build_rows(
        in_dirs=in_dirs,
        recursive=not args.no_recursive,
        sanitize=args.sanitize_sample_ids,
        sample_source=args.sample_source,
        jobs=args.jobs,
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
