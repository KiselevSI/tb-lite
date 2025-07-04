#!/usr/bin/env python3

import os
import argparse
import json
import re
import sys
import yaml
import gzip
from collections import defaultdict

class DrVariant:
    """Класс для хранения информации о варианте, вызывающем лекарственную устойчивость"""
    def __init__(self, gene_name, change, freq=1.0, drugs=None):
        self.gene_name = gene_name
        self.change = change
        self.freq = freq
        self.drugs = drugs or []
        self.pos = None
        self.effect = None
    
    def get_drugs(self):
        """Получить список препаратов, к которым вариант обеспечивает устойчивость"""
        return [d['drug'] for d in self.drugs]
    
    def __str__(self):
        return f"{self.gene_name} {self.change} ({self.freq:.2f})"

def parse_args():
    parser = argparse.ArgumentParser(description='Определение лекарственной устойчивости TB на основе аннотированных VCF файлов')
    parser.add_argument('-i', '--input', required=True, help='Путь к аннотированному VCF файлу')
    parser.add_argument('-o', '--output', required=True, help='Путь к выходному CSV файлу с отчетом (без расширения)')
    parser.add_argument('-d', '--detailed', action='store_true', help='Создать подробный отчет в формате TXT')
    parser.add_argument('-v', '--verbose', action='store_true', help='Подробный вывод')
    return parser.parse_args()

def load_gene_mapping(bed_path):
    """Загрузить отображение ID гена в имя и связи ген-препарат из BED файла"""
    gene_id2name = {}
    gene2drugs = {}
    
    with open(bed_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            row = line.strip().split('\t')
            if len(row) < 6:
                continue
                
            gene_id = row[3]  # ID гена в колонке 4
            gene_name = row[4]  # Имя гена в колонке 5
            drugs = row[5].split(',')  # Препараты в колонке 6
            
            gene_id2name[gene_id] = gene_name
            gene2drugs[gene_name] = drugs
    
    return gene_id2name, gene2drugs

def load_rules(rules_path):
    """Загрузить правила эпистаза из YAML файла"""
    if not rules_path or not os.path.exists(rules_path):
        return {}
    
    with open(rules_path, 'r') as f:
        rules = yaml.safe_load(f)
    
    return rules

def open_file(filename):
    """Открыть файл (сжатый или обычный)"""
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    return open(filename, 'r')

def parse_snpeff_annotation(ann_field):
    """Разобрать поле аннотации snpEff"""
    fields = ann_field.split('|')
    if len(fields) < 10:
        return None
    
    return {
        'allele': fields[0],
        'effect': fields[1],
        'impact': fields[2],
        'gene_id': fields[3],
        'gene_name': fields[4],
        'feature_type': fields[5],
        'feature_id': fields[6],
        'transcript_biotype': fields[7],
        'rank': fields[8],
        'hgvs_c': fields[9],
        'hgvs_p': fields[10] if len(fields) > 10 else "",
        'cdna_pos': fields[11] if len(fields) > 11 else "",
        'cds_pos': fields[12] if len(fields) > 12 else "",
        'protein_pos': fields[13] if len(fields) > 13 else "",
        'distance': fields[14] if len(fields) > 14 else "",
        'errors': fields[15] if len(fields) > 15 else ""
    }

def normalize_change(change):
    """Нормализовать представление мутации для лучшего сопоставления"""
    # Удалить общие префиксы
    change = re.sub(r'^(p\.|c\.|n\.)', '', change)
    
    # Обработать различные представления аминокислот
    aa_map = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
        'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
        'Stop': '*', '*': '*', 'del': 'del', 'ins': 'ins', 'fs': 'fs',
        'dup': 'dup'
    }
    
    # Конвертировать трехбуквенные коды аминокислот в однобуквенные
    for three_letter, one_letter in aa_map.items():
        change = change.replace(three_letter, one_letter)
    
    return change

def parse_vcf(vcf_file, gene_id2name):
    """Разобрать аннотированный с помощью snpEff VCF файл для извлечения вариантов"""
    variants = []
    
    try:
        with open_file(vcf_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                # Разобрать строку VCF
                cols = line.strip().split('\t')
                if len(cols) < 8:
                    continue
                    
                chrom = cols[0]
                pos = int(cols[1])
                variant_id = cols[2]
                ref = cols[3]
                alt = cols[4]
                qual = cols[5]
                filter_status = cols[6]
                info = cols[7]
                
                # Извлечь FORMAT и SAMPLE поля, если они доступны
                format_col = cols[8] if len(cols) > 8 else ""
                sample_col = cols[9] if len(cols) > 9 else ""
                
                # Извлечь частоту аллеля
                freq = 1.0
                # Попытаться получить AF из поля INFO
                af_match = re.search(r'AF=([^;]+)', info)
                if af_match:
                    freq = float(af_match.group(1))
                # Если нет AF в INFO, попытаться получить из поля образца
                elif format_col and sample_col and 'AD' in format_col:
                    ad_idx = format_col.split(':').index('AD')
                    ad_values = sample_col.split(':')[ad_idx].split(',')
                    if len(ad_values) >= 2 and sum(map(int, ad_values)) > 0:
                        ref_depth = int(ad_values[0])
                        alt_depth = int(ad_values[1])
                        total_depth = ref_depth + alt_depth
                        if total_depth > 0:
                            freq = alt_depth / total_depth
                
                # Извлечь аннотацию snpEff
                ann_match = re.search(r'ANN=([^;]+)', info)
                if not ann_match:
                    continue
                
                # Обработать все аннотации
                annotations = ann_match.group(1).split(',')
                
                for ann in annotations:
                    ann_data = parse_snpeff_annotation(ann)
                    if not ann_data:
                        continue
                    
                    gene_id = ann_data['gene_id']
                    gene_name = gene_id2name.get(gene_id, ann_data['gene_name'])
                    effect = ann_data['effect']
                    
                    # Определить строку изменения
                    change = ""
                    if ann_data['hgvs_p'] and ann_data['hgvs_p'] != "":
                        change = f"p.{ann_data['hgvs_p']}"
                    elif ann_data['hgvs_c'] and ann_data['hgvs_c'] != "":
                        change = f"c.{ann_data['hgvs_c']}"
                    else:
                        change = f"{ref}{pos}{alt}"
                    
                    # Извлечь информацию о глубине покрытия
                    depth = None
                    if format_col and sample_col:
                        if 'DP' in format_col:
                            dp_idx = format_col.split(':').index('DP')
                            depth = sample_col.split(':')[dp_idx]
                        elif 'AD' in format_col:
                            ad_idx = format_col.split(':').index('AD')
                            ad_values = sample_col.split(':')[ad_idx].split(',')
                            if len(ad_values) >= 2:
                                depth = sum(map(int, ad_values))
                    
                    # Создать объект варианта
                    variant = {
                        'chrom': chrom,
                        'pos': pos,
                        'ref': ref,
                        'alt': alt,
                        'effect': effect,
                        'gene_id': gene_id,
                        'gene_name': gene_name,
                        'change': change,
                        'norm_change': normalize_change(change),
                        'freq': freq,
                        'depth': depth,
                        'filter': filter_status,
                        'annotation': []  # Добавляем поле для аннотаций
                    }
                    
                    variants.append(variant)
    
    except Exception as e:
        print(f"Ошибка при разборе VCF файла {vcf_file}: {e}", file=sys.stderr)
        raise
    
    return variants

def identify_resistance(variants, resistance_db, gene2drugs):
    """Идентифицировать варианты лекарственной устойчивости"""
    dr_variants = []
    
    for variant in variants:
        gene_name = variant['gene_name']
        gene_id = variant['gene_id']
        change = variant['change']
        norm_change = variant['norm_change']
        
        # Проверить, есть ли ген в базе данных устойчивости
        for gene in [gene_name, gene_id]:
            if gene not in resistance_db:
                continue
            
            # Проверить на точное совпадение
            for db_change in resistance_db[gene]:
                db_norm_change = normalize_change(db_change)
                
                # Проверить на совпадение (точное или нормализованное)
                if (change == db_change or 
                    norm_change == db_norm_change or
                    norm_change in db_norm_change or
                    db_norm_change in norm_change):
                    
                    # Получить информацию о лекарственной устойчивости
                    drugs = []
                    for ann in resistance_db[gene][db_change]['annotations']:
                        if ann['type'] == 'drug_resistance':
                            drugs.append({
                                'drug': ann['drug'],
                                'confidence': ann.get('confidence', ''),
                                'comment': ann.get('comment', '')
                            })
                    
                    if drugs:
                        # Создать объект DrVariant
                        dr_variant = DrVariant(
                            gene_name=gene_name,
                            change=change,
                            freq=variant['freq'],
                            drugs=drugs
                        )
                        
                        # Добавляем позицию и эффект из исходного варианта
                        dr_variant.pos = variant['pos']
                        dr_variant.effect = variant['effect']
                        dr_variant.gene_id = variant['gene_id']
                        dr_variant.depth = variant['depth']
                        
                        dr_variants.append(dr_variant)
                        break  # Нашли совпадение, не нужно проверять другие мутации
            
            # Если нашли совпадение, не нужно проверять другой идентификатор гена
            if dr_variants and dr_variants[-1].gene_name == gene_name:
                break
        
        # Добавляем аннотации для других вариантов
        if gene_name in gene2drugs:
            drug_annotations = []
            for drug in gene2drugs[gene_name]:
                drug_annotations.append({
                    'type': 'who_confidence',
                    'drug': drug,
                    'confidence': 'Uncertain significance' if variant['effect'] == 'missense_variant' else 'Not assoc w R'
                })
            variant['annotation'] = drug_annotations
    
    return dr_variants

def apply_rules(dr_variants, rules):
    """Применить правила эпистаза для модификации предсказаний устойчивости"""
    # Это упрощенная реализация
    if not rules:
        return dr_variants
    
    # Обработать каждое правило
    for rule_name, rule in rules.items():
        if rule['type'] != 'epistasisRule':
            continue
        
        # Найти варианты исходного гена
        source_gene = rule['source']['gene_name']
        source_type = rule['source']['type']
        source_variants = [v for v in dr_variants if v.gene_name == source_gene]
        
        # Найти варианты целевого гена
        target_gene = rule['target']['gene_name']
        target_drugs = rule['target']['drug'] if isinstance(rule['target']['drug'], list) else [rule['target']['drug']]
        target_variants = [v for v in dr_variants if v.gene_name == target_gene]
        
        # Проверить на эпистаз
        if source_variants and target_variants:
            # Найти наивысшую частоту исходных вариантов
            max_source_freq = max([v.freq for v in source_variants])
            
            # Применить правило на основе частот
            if max_source_freq * 100 >= rule['source_inactivation_freq_cutoff']:
                # Модифицировать целевые варианты, чтобы удалить затронутые препараты
                for target_var in target_variants:
                    target_var.drugs = [d for d in target_var.drugs if d['drug'] not in target_drugs]
    
    # Удалить варианты без препаратов
    return [v for v in dr_variants if v.drugs]

def get_drtypes(dr_variants):
    """Определить тип устойчивости на основе препаратов"""
    resistant_drugs = set()
    for var in dr_variants:
        resistant_drugs.update(var.get_drugs())
    
    # Определить наборы препаратов для классификации
    FLQ_set = set(["levofloxacin", "moxifloxacin", "ciprofloxacin", "ofloxacin"])
    groupA_set = set(["bedaquiline", "linezolid"])
    
    # Проверить на специфические паттерны устойчивости
    rif = "rifampicin" in resistant_drugs
    inh = "isoniazid" in resistant_drugs
    flq = len(FLQ_set.intersection(resistant_drugs)) > 0
    gpa = len(groupA_set.intersection(resistant_drugs)) > 0
    
    # Классифицировать тип устойчивости
    if len(resistant_drugs) == 0:
        drtype = "Sensitive"
    elif (rif and not inh) and not flq:
        drtype = "RR-TB"
    elif (inh and not rif):
        drtype = "HR-TB"
    elif (rif and inh) and not flq:
        drtype = "MDR-TB"
    elif rif and (flq and not gpa):
        drtype = "Pre-XDR-TB"
    elif rif and (flq and gpa):
        drtype = "XDR-TB"
    else:
        drtype = "Other"
    
    return drtype

def generate_summary_report(sample_id, dr_variants, all_drugs, output_file):
    """Сгенерировать краткий отчет об устойчивости в формате CSV"""
    
    # Организовать варианты по препаратам
    drug_variants = defaultdict(list)
    for var in dr_variants:
        for drug in var.get_drugs():
            drug_variants[drug].append(var)
    
    # Сгенерировать отчет в формате CSV
    with open(output_file, 'w') as f:
        f.write("Drug,Genotypic Resistance,Mechanisms\n")
        
        # Фиксированный список препаратов в порядке, который хочет пользователь
        standard_drugs = [
            "rifampicin", "isoniazid", "ethambutol", "pyrazinamide", 
            "moxifloxacin", "levofloxacin", "bedaquiline", "delamanid", 
            "pretomanid", "linezolid", "streptomycin", "amikacin", 
            "kanamycin", "capreomycin", "clofazimine", "ethionamide",
            "para-aminosalicylic_acid", "cycloserine"
        ]
        
        # Убедиться, что все препараты в списке
        for drug in all_drugs:
            if drug not in standard_drugs:
                standard_drugs.append(drug)
        
        for drug in standard_drugs:
            if drug not in all_drugs:
                continue
                
            variants = drug_variants.get(drug, [])
            resistance = "R" if variants else ""
            mechanisms = ", ".join([f"{v.gene_name} {v.change} ({v.freq:.2f})" for v in variants])
            
            drug_name = drug.replace('_', '-')
            f.write(f"{drug_name},{resistance},{mechanisms}\n")

def generate_detailed_report(sample_id, dr_variants, other_variants, all_drugs, output_file):
    """Сгенерировать подробный отчет об устойчивости в формате TXT"""
    
    # Организовать варианты по препаратам
    drug_variants = defaultdict(list)
    for var in dr_variants:
        for drug_info in var.drugs:
            drug = drug_info['drug']
            drug_variants[drug].append((var, drug_info))
    
    # Исправляем дублирование префиксов в именах мутаций
    for var in dr_variants + other_variants:
        if var.change.startswith('p.p.'):
            var.change = var.change.replace('p.p.', 'p.')
        elif var.change.startswith('c.c.'):
            var.change = var.change.replace('c.c.', 'c.')
        elif var.change.startswith('n.n.'):
            var.change = var.change.replace('n.n.', 'n.')
    
    # Формируем отчет
    with open(output_file, 'w') as f:
        # Первая секция - Resistance report (краткая сводка)
        f.write("Resistance report\n")
        f.write("-----------------\n")
        f.write("Drug\tGenotypic Resistance\tMechanisms\n")
        
        # Фиксированный список препаратов
        standard_drugs = [
            "rifampicin", "isoniazid", "ethambutol", "pyrazinamide", 
            "moxifloxacin", "levofloxacin", "bedaquiline", "delamanid", 
            "pretomanid", "linezolid", "streptomycin", "amikacin", 
            "kanamycin", "capreomycin", "clofazimine", "ethionamide",
            "para-aminosalicylic_acid", "cycloserine"
        ]
        
        # Убедиться, что все препараты в списке
        for drug in all_drugs:
            if drug not in standard_drugs:
                standard_drugs.append(drug)
        
        # Генерация первой секции с краткой сводкой
        drug_resistance_summary = {}
        for drug in standard_drugs:
            if drug not in all_drugs:
                continue
                
            variants_list = []
            for var in dr_variants:
                if drug in var.get_drugs():
                    variants_list.append(var)
            
            resistance = "R" if variants_list else ""
            mechanisms = ", ".join([f"{var.gene_name} {var.change} ({var.freq:.2f})" for v in variants_list])
            
            drug_name = drug.replace('_', ' ').title()
            f.write(f"{drug_name}\t{resistance}\t{mechanisms}\n")
            
            # Сохраняем для использования ниже
            drug_resistance_summary[drug] = (resistance, mechanisms)
        
        # Вторая секция - Resistance variants report
        f.write("\nResistance variants report\n")
        f.write("-----------------\n")
        f.write("Genome Position\tGene Name\tLocus Tag\tVariant Type\tEstimated Fraction\tDrug\tConfidence\tComment\n")
        
        # Сортируем варианты по позиции в геноме
        sorted_variants = sorted(dr_variants, key=lambda x: x.pos if hasattr(x, 'pos') else 0)
        
        for var in sorted_variants:
            # Собираем информацию о препаратах для этого варианта
            drugs = []
            confidences = []
            comments = []
            
            for drug_info in var.drugs:
                drugs.append(drug_info['drug'])
                confidences.append(drug_info.get('confidence', ''))
                if 'comment' in drug_info and drug_info['comment']:
                    comments.append(drug_info['comment'])
            
            # Формируем строки для отчета
            drugs_str = ",".join(drugs)
            confidences_str = ",".join(filter(None, confidences))
            comments_str = ",".join(filter(None, comments))
            
            # Получаем дополнительную информацию
            pos = getattr(var, 'pos', '')
            locus_tag = getattr(var, 'gene_id', '')
            variant_type = getattr(var, 'effect', '')
            
            f.write(f"{pos}\t{var.gene_name}\t{locus_tag}\t{variant_type}\t{var.freq:.3f}\t{drugs_str}\t{confidences_str}\t{comments_str}\n")
        
        # Третья секция - Other variants report
        f.write("\nOther variants report\n")
        f.write("---------------------\n")
        f.write("Genome Position\tGene Name\tLocus Tag\tVariant Type\tEstimated Fraction\tGene Associated Drug\tConfidence\tComment\n")
        
        # Сортируем другие варианты по позиции в геноме
        sorted_other_variants = sorted(other_variants, key=lambda x: x.pos if hasattr(x, 'pos') else 0)
        
        for var in sorted_other_variants:
            # Проверяем наличие аннотаций
            annotations = getattr(var, 'annotation', [])
            if not annotations or not any(a.get('drug') for a in annotations):
                continue
            
            # Получаем основную информацию о варианте
            pos = getattr(var, 'pos', '')
            locus_tag = getattr(var, 'gene_id', '')
            variant_type = getattr(var, 'effect', '')
            
            # Обрабатываем аннотации
            first_row = True
            for ann in annotations:
                if 'drug' not in ann:
                    continue
                
                drug = ann['drug']
                confidence = ann.get('confidence', '')
                comment = ann.get('comment', '')
                f.write(f"{pos}\t{var.gene_name}\t{locus_tag}\t{variant_type}\t{var.freq:.3f}\t{drug}\t{confidence}\t{comment}\n")
                    
                #if first_row:
                #    f.write(f"{pos}\t{var.gene_name}\t{locus_tag}\t{variant_type}\t{var.freq:.3f}\t{drug}\t{confidence}\t{comment}\n")
                #    first_row = False
                #else:
                #    f.write(f"\t\t\t\t\t{drug}\t{confidence}\t{comment}\n")

def main():
    args = parse_args()
    
    # Проверяем существование входного файла
    if not os.path.exists(args.input):
        print(f"Ошибка: Входной файл {args.input} не найден", file=sys.stderr)
        sys.exit(1)
    
    # Проверяем, что входной файл - это действительно файл, а не директория
    if os.path.isdir(args.input):
        print(f"Ошибка: {args.input} является директорией, а не файлом. Укажите путь к VCF файлу.", file=sys.stderr)
        sys.exit(1)
    
    # Проверка расширения файла
    if not (args.input.endswith('.vcf') or args.input.endswith('.vcf.gz')):
        print(f"Предупреждение: Файл {args.input} не имеет расширения .vcf или .vcf.gz", file=sys.stderr)
    
    # Определяем пути для выходных файлов
    csv_output = f"{args.output}"
    txt_output = args.output.replace('.csv', '.detailed.txt') if args.detailed else None
    
    # Создаем директорию для выходных файлов, если нужно
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # Определить пути к файлам базы данных
    script_dir = os.path.dirname(os.path.abspath(__file__))
    db_path = os.path.join(script_dir, 'db_drug_resist', 'tbdb.dr.json')
    bed_path = os.path.join(script_dir, 'db_drug_resist', 'tbdb.bed')
    rules_path = os.path.join(script_dir, 'db_drug_resist', 'tbdb.rules.yml')
    
    # Загрузить файлы базы данных
    try:
        with open(db_path, 'r') as f:
            resistance_db = json.load(f)
        
        gene_id2name, gene2drugs = load_gene_mapping(bed_path)
        rules = load_rules(rules_path) if os.path.exists(rules_path) else {}
        
        # Получить список всех препаратов
        all_drugs = set()
        for drugs in gene2drugs.values():
            all_drugs.update(drugs)
        
        if args.verbose:
            print(f"Загружена база данных устойчивости с {len(resistance_db)} генами")
            print(f"Загружено отображение генов с {len(gene_id2name)} генами")
            print(f"Найдено {len(all_drugs)} препаратов в базе данных")
            if rules:
                print(f"Загружено {len(rules)} правил устойчивости")
    
    except FileNotFoundError as e:
        print(f"Ошибка: Файл базы данных не найден: {e}", file=sys.stderr)
        print(f"Проверьте, что в директории {os.path.join(script_dir, 'db')} есть необходимые файлы:", file=sys.stderr)
        print(f"- tbdb.dr.json: база данных мутаций устойчивости", file=sys.stderr)
        print(f"- tbdb.bed: файл сопоставления генов и препаратов", file=sys.stderr)
        print(f"- tbdb.rules.yml: (опционально) правила эпистаза", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Ошибка при загрузке файлов базы данных: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Получить имя образца из имени файла
    sample_id = os.path.basename(args.input).split('.')[0]
    
    try:
        if args.verbose:
            print(f"Обработка файла {args.input}...")
        
        # Анализ VCF файла
        variants = parse_vcf(args.input, gene_id2name)
        if args.verbose:
            print(f"Найдено {len(variants)} вариантов")
        
        # Идентификация вариантов устойчивости
        dr_variants = identify_resistance(variants, resistance_db, gene2drugs)
        if rules:
            dr_variants = apply_rules(dr_variants, rules)
        
        if args.verbose:
            print(f"Найдено {len(dr_variants)} вариантов устойчивости")
        
        # Выделение других вариантов (не связанных с устойчивостью)
        other_variants = []
        dr_genes_positions = {(v.gene_name, v.change) for v in dr_variants}
        for v in variants:
            if (v['gene_name'], v['change']) not in dr_genes_positions:
                # Создаем объект для других вариантов с нужными полями
                var_obj = type('NonDrVariant', (), {
                    'gene_name': v['gene_name'],
                    'change': v['change'],
                    'pos': v.get('pos', ''),
                    'gene_id': v.get('gene_id', ''),
                    'effect': v.get('effect', ''),
                    'depth': v.get('depth', ''),
                    'freq': v.get('freq', 0.0),
                    'annotation': v.get('annotation', [])
                })
                other_variants.append(var_obj)
        
        # Генерация отчетов
        generate_summary_report(sample_id, dr_variants, all_drugs, csv_output)
        if args.detailed:
            generate_detailed_report(sample_id, dr_variants, other_variants, all_drugs, txt_output)
        
        if args.verbose:
            print(f"Краткий отчет об устойчивости сгенерирован: {csv_output}")
            if args.detailed:
                print(f"Подробный отчет об устойчивости сгенерирован: {txt_output}")
    
    except Exception as e:
        print(f"Ошибка при обработке {args.input}: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()