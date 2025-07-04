#!/usr/bin/env python

import os
import argparse
import json
import glob
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
    
    def get_drugs(self):
        """Получить список препаратов, к которым вариант обеспечивает устойчивость"""
        return [d['drug'] for d in self.drugs]
    
    def __str__(self):
        return f"{self.gene_name} {self.change} ({self.freq:.2f})"

def parse_args():
    parser = argparse.ArgumentParser(description='Определение лекарственной устойчивости TB на основе аннотированных VCF файлов')
    parser.add_argument('-i', '--input', required=True, help='Директория с аннотированными VCF файлами')
    parser.add_argument('-o', '--output', required=True, help='Директория для выходных отчетов')
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
                        change = ann_data['hgvs_p']
                    elif ann_data['hgvs_c'] and ann_data['hgvs_c'] != "":
                        change = ann_data['hgvs_c']
                    else:
                        change = f"{ref}{pos}{alt}"
                    
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
                        'filter': filter_status
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
                        
                        dr_variants.append(dr_variant)
                        break  # Нашли совпадение, не нужно проверять другие мутации
            
            # Если нашли совпадение, не нужно проверять другой идентификатор гена
            if dr_variants and dr_variants[-1].gene_name == gene_name:
                break
    
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

def generate_report(sample_id, dr_variants, all_drugs, output_dir):
    """Сгенерировать отчет об устойчивости в формате TSV"""
    report_file = os.path.join(output_dir, f"{sample_id}.resistance.tsv")
    
    # Организовать варианты по препаратам
    drug_variants = defaultdict(list)
    for var in dr_variants:
        for drug in var.get_drugs():
            drug_variants[drug].append(var)
    
    # Сгенерировать отчет
    with open(report_file, 'w') as f:
        f.write("Drug\tResistance\tMechanisms\n")
        
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
            
            drug_name = drug.replace('_', ' ').title()
            f.write(f"{drug_name}\t{resistance}\t{mechanisms}\n")
        
        # Добавить тип устойчивости
        drtype = get_drtypes(dr_variants)
        f.write(f"\nResistance Type\t{drtype}\t\n")

def main():
    args = parse_args()
    
    # Убедиться, что выходная директория существует
    os.makedirs(args.output, exist_ok=True)
    
    # Определить пути к файлам базы данных
    script_dir = os.path.dirname(os.path.abspath(__file__))
    db_path = os.path.join(script_dir, 'db', 'tbdb.dr.json')
    bed_path = os.path.join(script_dir, 'db', 'tbdb.bed')
    rules_path = os.path.join(script_dir, 'db', 'tbdb.rules.yml')
    
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
    
    except Exception as e:
        print(f"Ошибка при загрузке файлов базы данных: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Обработать каждый VCF файл
    vcf_files = glob.glob(os.path.join(args.input, '*.vcf*'))
    if not vcf_files:
        print(f"VCF файлы не найдены в {args.input}", file=sys.stderr)
        sys.exit(1)
    
    for vcf_file in vcf_files:
        sample_id = os.path.basename(vcf_file).split('.')[0]
        if args.verbose:
            print(f"Обработка {sample_id}...")
        
        try:
            variants = parse_vcf(vcf_file, gene_id2name)
            if args.verbose:
                print(f"Найдено {len(variants)} вариантов в {sample_id}")
            
            dr_variants = identify_resistance(variants, resistance_db, gene2drugs)
            if rules:
                dr_variants = apply_rules(dr_variants, rules)
            
            if args.verbose:
                print(f"Найдено {len(dr_variants)} вариантов устойчивости в {sample_id}")
            
            generate_report(sample_id, dr_variants, all_drugs, args.output)
            
            if args.verbose:
                print(f"Отчет об устойчивости сгенерирован для {sample_id}")
        
        except Exception as e:
            print(f"Ошибка при обработке {vcf_file}: {e}", file=sys.stderr)
            continue

if __name__ == "__main__":
    main()