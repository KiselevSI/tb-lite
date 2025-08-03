#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
from collections import defaultdict
import tempfile
import shutil

import pysam


def merge_vcfs_simple(vcf_paths, out_path):
    """
    Слияние VCF по ключу (CHROM, POS, REF):
      - ALT = объединение всех альтернатив без дубликатов (порядок стабилизирован)
      - QUAL = максимум
      - Остальные поля берутся из записи с максимальным QUAL

    Важно: генотипы/FORMAT и INFO НЕ сводятся между файлами.
    """
    if not vcf_paths:
        raise ValueError("Не переданы входные файлы VCF.")

    # Шаблон заголовка берём из первого файла
    vcf_in0 = pysam.VariantFile(vcf_paths[0])
    header = vcf_in0.header.copy()
    vcf_in0.close()

    # Добавим мета-строку о происхождении файла
    # (pysam поддерживает add_meta / add_line для модификации заголовка)
    header.add_meta(
        "source",
        items=[("Value", "vcf-merge-simple"), ("Description", "Union ALT, max QUAL per (CHROM,POS,REF)")]
    )

    # Контроль порядка контигов для сортировки вывода
    contig_order = {name: i for i, name in enumerate(header.contigs.keys())}

    # Копим кандидатов по ключу
    buckets = defaultdict(list)

    for path in vcf_paths:
        vcf_in = pysam.VariantFile(path)
        for rec in vcf_in.fetch():
            key = (rec.contig, rec.pos, rec.ref)
            # храним ссылку на запись и подготовленные значения
            buckets[key].append(rec)
        vcf_in.close()

    # Пишем во временный (не сжатый) файл; потом при необходимости bgzip+tabix
    write_to_gz = out_path.endswith(".vcf.gz")
    if write_to_gz:
        tmp_vcf = tempfile.NamedTemporaryFile(prefix="merge_tmp_", suffix=".vcf", delete=False).name
        dst_vcf_path = tmp_vcf
    else:
        dst_vcf_path = out_path

    vcf_out = pysam.VariantFile(dst_vcf_path, "w", header=header)

    # Стабильная сортировка по порядку контигов из заголовка, затем по POS
    def contig_rank(chrom):
        # неизвестные контиги отправляем в конец
        return contig_order.get(chrom, len(contig_order))

    for key in sorted(buckets.keys(), key=lambda k: (contig_rank(k[0]), k[1], k[2])):
        records = buckets[key]
        if len(records) == 1:
            vcf_out.write(records[0])
            continue

        # Выбираем запись с максимальным QUAL как «донор»
        best = max(records, key=lambda r: (r.qual if r.qual is not None else float("-inf")))

        # Собираем ALT из всех записей
        alt_set = []
        seen = set()
        for r in records:
            if r.alts:
                for a in r.alts:
                    if a not in seen:
                        seen.add(a)
                        alt_set.append(a)

        # Создаём новый VariantRecord на базе best
        new_rec = best.copy()
        # alleles = (REF, ALT1, ALT2, ...)
        new_rec.alleles = tuple([new_rec.ref] + alt_set)
        # QUAL = максимум
        best_qual = max((r.qual for r in records if r.qual is not None), default=None)
        new_rec.qual = best_qual

        # Важно: FORMAT/генотипы и INFO остаются как у best.
        # Если нужно объединять INFO — добавьте свою логику здесь.

        vcf_out.write(new_rec)

    vcf_out.close()

    # Если пользователь ждёт .vcf.gz — bgzip + tabix, как рекомендует pysam
    if write_to_gz:
        # bgzip + tabix: pysam.tabix_index умеет сам сжимать и индексировать текстовый VCF
        # возвращает путь к .gz
        gz_path = pysam.tabix_index(
            dst_vcf_path, preset="vcf", force=True, keep_original=False
        )
        # Переместим в требуемое имя
        shutil.move(gz_path, out_path)
        # Индекс уже создан рядом (.tbi)
    return out_path


def main():
    parser = argparse.ArgumentParser(
        description="Объединение VCF по (CHROM, POS, REF) с объединением ALT и максимальным QUAL."
    )
    parser.add_argument(
        "-v", "--vcf", nargs="+", required=True,
        help="Список входных VCF (.vcf или .vcf.gz). Порядок не важен."
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Выходной файл (.vcf или .vcf.gz). Для .vcf.gz будет создан индекс .tbi."
    )
    args = parser.parse_args()

    try:
        merge_vcfs_simple(args.vcf, args.output)
    except Exception as e:
        print(f"Ошибка: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
