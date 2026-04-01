#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path

# шаблоны для поиска номеров ридов
R1_PATTERNS = [r'(_R?1)(?=[._])', r'(_1)(?=[._])']
R2_PATTERNS = [r'(_R?2)(?=[._])', r'(_2)(?=[._])']

FASTQ_EXTS = ('.fastq.gz', '.fq.gz')
PLAIN_FASTQ_EXTS = ('.fastq', '.fq')

def detect_read_label(name: str):
    """Вернёт ('R1'|'R2'|None, core_name) где core_name — база для sample id."""
    for p1, p2 in zip(R1_PATTERNS, R2_PATTERNS):
        m1 = re.search(p1, name, flags=re.IGNORECASE)
        m2 = re.search(p2, name, flags=re.IGNORECASE)
        if m1:
            core = name[:m1.start()] + name[m1.end():]
            return 'R1', core
        if m2:
            core = name[:m2.start()] + name[m2.end():]
            return 'R2', core
    return None, name  # одиночные

def sample_id_from_core(core: str):
    """Убираем расширения и лишние точки/подчёркивания на конце."""
    core = re.sub(r'\.(fastq|fq)(\.gz)?$', '', core, flags=re.IGNORECASE)
    return core.rstrip('._-')

def build_samplesheet(in_dir: Path):
    files = [p for p in in_dir.rglob('*') if p.is_file() and p.suffixes and ''.join(p.suffixes).lower() in FASTQ_EXTS]
    plain_fastq = [p for p in in_dir.rglob('*') if p.is_file() and p.suffixes and ''.join(p.suffixes).lower() in PLAIN_FASTQ_EXTS]

    if not files and plain_fastq:
        raise SystemExit(
            "ERROR: найдены только несжатые FASTQ/FQ файлы. "
            "TB-Lite ожидает gzipped входы (*.fastq.gz или *.fq.gz). "
            "Сначала сожмите файлы, затем создайте samplesheet."
        )
    if files and plain_fastq:
        print(
            f"WARNING: найдено {len(plain_fastq)} несжатых FASTQ/FQ файлов. "
            "Они будут проигнорированы; в samplesheet попадут только *.fastq.gz/*.fq.gz."
        )
    pairs = {}
    singles = []

    for f in sorted(files):
        label, core = detect_read_label(f.name)
        sid = sample_id_from_core(core)
        if label == 'R1':
            pairs.setdefault(sid, {})['R1'] = str(f)
        elif label == 'R2':
            pairs.setdefault(sid, {})['R2'] = str(f)
        else:
            singles.append((sid, str(f)))

    rows = []

    # сначала пары
    for sid, info in pairs.items():
        r1 = info.get('R1')
        r2 = info.get('R2')
        if r1 and r2:
            rows.append([sid, r1, r2])
        elif r1 and not r2:
            rows.append([sid, r1, ''])
        elif r2 and not r1:
            rows.append([sid, '', r2])  # экзотичный случай

    # одиночные, если их sample_id ещё не встречался в парах
    for sid, r1 in singles:
        if sid not in pairs:  # иначе уже добавили
            rows.append([sid, r1, ''])

    # сортируем по sample id для стабильности
    rows.sort(key=lambda x: x[0])
    return rows

def main():
    ap = argparse.ArgumentParser(description='Создать samplesheet.csv для Nextflow из директории с gzipped FASTQ.')
    ap.add_argument('-i', '--input', required=True, help='Входная директория с файлами *.fastq.gz или *.fq.gz')
    ap.add_argument('-o', '--output', required=True, help='Путь к выходному CSV')
    args = ap.parse_args()

    in_dir = Path(args.input).resolve()
    if not in_dir.is_dir():
        raise SystemExit(f"ERROR: '{in_dir}' не директория")

    rows = build_samplesheet(in_dir)

    out_path = Path(args.output).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open('w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['sample', 'fastq_1', 'fastq_2'])
        w.writerows(rows)

    print(f"✓ Готово: {out_path} (samples: {len(rows)})")

if __name__ == '__main__':
    main()
