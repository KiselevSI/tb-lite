#!/usr/bin/env python3
"""
pileup_counter.py  -- глубина **после фильтров** + фильтрованные средние BQ/MAPQ
"""
import pysam, argparse, gzip, csv
from collections import Counter

def read_regions(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 3:
                continue
            yield row[0], int(row[1]), int(row[2])

def count_site(bam, fasta, chrom, pos, mapq_thr, bq_thr):
    depth_raw  = 0      # общая глубина (оставим, вдруг понадобится)
    depth_filt = 0      # глубина после фильтров – она пойдёт в отчёт
    qual_sum   = 0
    mapq_sum   = 0
    counts     = Counter()

    for col in bam.pileup(chrom, pos-1, pos,
                          min_base_quality=0,
                          min_mapping_quality=0,
                          truncate=True,
                          stepper="nofilter"):
        if col.reference_pos != pos-1:
            continue
        for read in col.pileups:
            if read.is_del or read.is_refskip:
                continue

            depth_raw += 1                     # считаем для статистики

            base   = read.alignment.query_sequence[read.query_position]
            strand = "-" if read.alignment.is_reverse else "+"
            bq_val = read.alignment.query_qualities[read.query_position]
            mq_val = read.alignment.mapping_quality

            if mq_val >= mapq_thr and bq_val >= bq_thr:
                counts[(base, strand)] += 1
                depth_filt += 1
                qual_sum   += bq_val
                mapq_sum   += mq_val

    ref_base = fasta.fetch(chrom, pos-1, pos).upper()
    avg_bq   = round(qual_sum / depth_filt, 2) if depth_filt else 0
    avg_mq   = round(mapq_sum / depth_filt, 2) if depth_filt else 0
    return ref_base, depth_filt, counts, avg_bq, avg_mq   # <-- depth_filt!

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-b", "--bam",   required=True)
    ap.add_argument("-f", "--fasta", required=True)
    ap.add_argument("-l", "--list",  required=True)
    ap.add_argument("-q", "--mapq",  type=int, default=20)
    ap.add_argument("-Q", "--bq",    type=int, default=20)
    args = ap.parse_args()

    bam   = pysam.AlignmentFile(args.bam, "rb")
    fasta = pysam.FastaFile(args.fasta)

    header = ["chrom","pos","ref","depth",   # заголовок поменяли
              "avg_bq","avg_mapq",
              "A","C","G","T","N"]
    print("\t".join(header))

    for chrom, start, end in read_regions(args.list):
        for pos in range(start, end+1):
            ref, depth, c, avg_bq, avg_mq = count_site(
                bam, fasta, chrom, pos, args.mapq, args.bq
            )
            row = [chrom, pos, ref, depth, avg_bq, avg_mq]
            for base in "ACGTN":
                row.append(c.get((base, "+"), 0) + c.get((base, "-"), 0))
            print("\t".join(map(str, row)))

if __name__ == "__main__":
    main()
