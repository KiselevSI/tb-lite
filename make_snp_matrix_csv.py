#!/usr/bin/env python3
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import gzip
import os
import re
import sys
import tempfile
from pathlib import Path

VCF_EXTS = ('.vcf', '.vcf.gz')
VALID_SAMPLE_RE = re.compile(r'^[A-Za-z0-9._-]+$')
PROGRESS_INTERVAL = 1000


def is_vcf(path: Path) -> bool:
    name = path.name.lower()
    return any(name.endswith(ext) for ext in VCF_EXTS)


def is_vcf_name(name: str) -> bool:
    name = name.lower()
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
            files.append(path.absolute())
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


def iter_vcfs(root: Path, recursive: bool):
    if root.is_file():
        if is_vcf_name(root.name):
            yield root.absolute()
        return

    if not root.is_dir():
        return

    stack = [str(root)]
    while stack:
        current = stack.pop()
        with os.scandir(current) as entries:
            for entry in entries:
                try:
                    if recursive and entry.is_dir(follow_symlinks=False):
                        stack.append(entry.path)
                    elif entry.is_file(follow_symlinks=True) and is_vcf_name(entry.name):
                        yield Path(entry.path).absolute()
                except OSError as exc:
                    raise ValueError(f"ERROR: cannot read '{entry.path}': {exc}") from exc


def prepare_scan_units(in_dirs: list[Path], recursive: bool, jobs: int) -> list[list[Path]]:
    units = []
    target_jobs = max(1, jobs)

    for in_dir in in_dirs:
        entries = []
        for count, entry in enumerate(in_dir.iterdir(), start=1):
            entries.append(entry)
            if count % PROGRESS_INTERVAL == 0:
                print(f"Preparing scan {in_dir}: seen {count} entries", file=sys.stderr)
        print(f"Preparing scan {in_dir}: seen {len(entries)} entries", file=sys.stderr)

        if not entries:
            units.append([in_dir])
            continue

        dirs = [entry for entry in entries if recursive and entry.is_dir()]
        files = [entry for entry in entries if entry.is_file()]

        units.extend([directory] for directory in dirs)

        if files:
            chunk_size = max(1, (len(files) + target_jobs - 1) // target_jobs)
            units.extend(files[index:index + chunk_size] for index in range(0, len(files), chunk_size))

        if not recursive:
            other_entries = [entry for entry in entries if not entry.is_file()]
            if not files and other_entries:
                units.append([in_dir])

    return units


def distribute_units(units: list[list[Path]], jobs: int) -> list[list[Path]]:
    worker_count = min(max(1, jobs), max(1, len(units)))
    shards = [[] for _ in range(worker_count)]
    for index, unit in enumerate(units):
        shards[index % worker_count].extend(unit)
    return [shard for shard in shards if shard]


def write_filename_shard(shard_id: int, roots: list[Path], recursive: bool, sanitize: bool, tmp_dir: Path) -> tuple[Path, dict[str, str], int]:
    shard_path = tmp_dir / f"part-{shard_id:04d}.csv"
    seen_samples = {}
    count = 0

    with shard_path.open('w', newline='') as fh:
        writer = csv.writer(fh)
        for root in roots:
            for vcf in iter_vcfs(root, recursive):
                sample = sample_id_from_vcf(vcf)
                if sanitize:
                    sample = sanitize_sample_id(sample)
                if not sample:
                    raise ValueError(f"ERROR: cannot derive sample ID from '{vcf}'")
                if not VALID_SAMPLE_RE.fullmatch(sample):
                    raise ValueError(
                        f"ERROR: sample ID '{sample}' from '{vcf.name}' contains unsupported characters. "
                        "Rename the file or rerun with --sanitize-sample-ids."
                    )
                vcf_str = str(vcf)
                if sample in seen_samples:
                    raise ValueError(
                        f"ERROR: duplicate sample ID '{sample}' derived from both "
                        f"'{seen_samples[sample]}' and '{vcf_str}'"
                    )
                seen_samples[sample] = vcf_str
                writer.writerow([sample, vcf_str])
                count += 1
                if count % PROGRESS_INTERVAL == 0:
                    print(f"Worker {shard_id}: wrote {count} VCF rows", file=sys.stderr)

    print(f"Worker {shard_id}: wrote {count} VCF rows", file=sys.stderr)
    return shard_path, seen_samples, count


def merge_filename_shards(shards: list[tuple[Path, dict[str, str], int]], out_path: Path) -> int:
    seen_samples = {}
    total = 0
    expected_total = sum(item[2] for item in shards)

    for _shard_path, shard_seen, _count in shards:
        for sample, vcf in shard_seen.items():
            if sample in seen_samples:
                raise SystemExit(
                    f"ERROR: duplicate sample ID '{sample}' derived from both "
                    f"'{seen_samples[sample]}' and '{vcf}'"
                )
            seen_samples[sample] = vcf

    tmp_output = out_path.with_name(f".{out_path.name}.tmp")
    with tmp_output.open('w', newline='') as out_fh:
        writer = csv.writer(out_fh)
        writer.writerow(['sample', 'vcf'])
        for shard_path, _shard_seen, count in shards:
            with shard_path.open(newline='') as shard_fh:
                for row in csv.reader(shard_fh):
                    writer.writerow(row)
            total += count
            report_progress(total, expected_total)

    tmp_output.replace(out_path)
    return total


def write_filename_csv(in_dirs: list[Path], recursive: bool, sanitize: bool, jobs: int, out_path: Path) -> int:
    units = prepare_scan_units(in_dirs, recursive, jobs)
    shards = distribute_units(units, jobs)

    if not shards:
        raise SystemExit("ERROR: VCF files (*.vcf or *.vcf.gz) not found in input directories")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{out_path.name}.parts.", dir=out_path.parent) as tmp_name:
        tmp_dir = Path(tmp_name)
        shard_results = []
        with ThreadPoolExecutor(max_workers=len(shards)) as executor:
            futures = [
                executor.submit(write_filename_shard, shard_id, shard_roots, recursive, sanitize, tmp_dir)
                for shard_id, shard_roots in enumerate(shards, start=1)
            ]
            for future in as_completed(futures):
                try:
                    shard_results.append(future.result())
                except ValueError as exc:
                    raise SystemExit(str(exc)) from exc

        if not any(count for _path, _seen, count in shard_results):
            raise SystemExit("ERROR: VCF files (*.vcf or *.vcf.gz) not found in input directories")

        shard_results.sort(key=lambda item: item[0].name)
        return merge_filename_shards(shard_results, out_path)


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
        help='Number of parallel workers. Header mode reads VCF headers in parallel; filename mode scans/writes shard files in parallel. Default: 1',
    )
    args = ap.parse_args()

    if args.jobs < 1:
        raise SystemExit("ERROR: --jobs must be a positive integer")

    in_dirs = [Path(path).resolve() for path in args.input]
    for in_dir in in_dirs:
        if not in_dir.is_dir():
            raise SystemExit(f"ERROR: '{in_dir}' is not a directory")

    out_path = Path(args.output).absolute()
    if args.sample_source == 'filename':
        total = write_filename_csv(
            in_dirs=in_dirs,
            recursive=not args.no_recursive,
            sanitize=args.sanitize_sample_ids,
            jobs=args.jobs,
            out_path=out_path,
        )
        print(f"Done: {out_path} (samples: {total})")
        return

    rows = build_rows(
        in_dirs=in_dirs,
        recursive=not args.no_recursive,
        sanitize=args.sanitize_sample_ids,
        sample_source=args.sample_source,
        jobs=args.jobs,
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open('w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['sample', 'vcf'])
        writer.writerows(rows)

    print(f"Done: {out_path} (samples: {len(rows)})")


if __name__ == '__main__':
    main()
