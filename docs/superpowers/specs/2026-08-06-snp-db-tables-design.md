# SNP-таблицы для TB Platform из пайплайна tb-lite

Дата: 2026-08-06
Статус: утверждён, готов к реализации

## Задача

После прогона tb-lite на 200 000+ образцах нужно уметь добавить их на сайт TB Platform.
Сейчас пайплайн отдаёт таблицы для `general`, `general_lineage`, `general_spoligo`,
`general_spoligo_spacers`, `deletions`, `general_drug_resist`, `tb_mix_lineage` и
`is6110_*`, но **не отдаёт ничего для VCF**. Без этого не работают:

- экспорт SNP matrix TSV (`snp_sites` + `sample_snp_alleles`);
- поиск по SNP и smart search (`snp_sites` + `sample_snp_alleles`, `vcf_table`);
- feature analysis / clade-признаки (`vcf_table`);
- построение дерева из БД и grafting (`sample_snp_profiles`).

Пайплайн должен сгенерировать файлы, готовые к импорту. Импорт запускается **вручную**,
пайплайн к БД не подключается.

## Проверенные факты

Эти факты установлены измерением, а не предположением, и на них опирается весь дизайн.

### site_key воспроизводит прод

`site_key` в проде — sha1-hex от `"\x1f".join([chrom, pos, ref, alt, effect, impact, gene,
gene_id, product_name, strand, feature_id, biotype, rank, hgvs_c, hgvs_p, cds_pos, aa_pos])`.
Пересчёт по этой формуле для трёх строк из `tb_database.dump` даёт побитовое совпадение:

| pos | site_key (прод и пересчёт) |
|---|---|
| 1977 | `7229c7508aa7fb285706a91c398ca4b22bdccee9` |
| 4013 | `f03d041b126f9db3c71d4f62fb698395201f8ee5` |
| 8342 | `a9aa271b2a7b3c733c1a9b2cd77fbfa3911ca974` |

Следствие: новые данные корректно схлопываются с существующими 4.38M строк `snp_sites`
через `ON CONFLICT (site_key)`. Любое изменение набора полей или их нормализации ломает
совместимость — формулу трогать нельзя.

### Семантика значений в проде

- `snp_sites.product_name` и `snp_sites.strand` — **NULL**. В per-sample VCF нет INFO-полей
  `name`/`strand`, они появляются только в `FINAL_ANNOTATION_TABLE.tsv` из
  `POST_PROCESS_TABLE` по feature table. Заполнять их нельзя — изменится `site_key`.
- `snp_sites.qual` — максимум QUAL по всем образцам, у которых встретился этот `site_key`.
- `sample_snp_alleles.filter` — NULL (в VCF `FILTER = "."`).
- `sample_snp_alleles.af` — 1 (INFO/AF от freebayes при `--ploidy 1`).
- `sample_snp_alleles.dp` — FORMAT/DP, `qual` — QUAL записи.
- `vcf_table` содержит indel'ы (пример из прода: `SRR9738528 1552 TAAAAAA`).

### Пишем всё содержимое VCF

Никакой фильтрации при извлечении: ни `--snps-only`, ни `--pass-only`, ни `--min-dp`,
ни `--min-af`, ни масок. Проверка на реальном файле прогона (ERR4797591, 1589 записей):
1423 `TYPE=snp`, 64 `del`, 51 `complex`, 48 `ins`, 3 `mnp`, все с `GT=1`.

VCF уже жёстко отфильтрован на этапе вызова вариантов, повторно фильтровать нечего:
`freebayes --ploidy 1 --min-mapping-quality 20 --min-base-quality 20 --min-coverage 10
--min-alternate-fraction 0.80`, затем `BCFTOOLS_VIEW -i 'QUAL>=10 && INFO/DP>=10'`.

### Источник данных

`results/annotate_vcf/<sample>/<sample>.annotated.ann.vcf` — per-sample VCF после snpEff
(`SNPEFF_SNPEFF` в `CALLVAR`), ровно одна колонка образца, `FORMAT = GT:DP:AD:RO:QR:AO:QA:GL`.
В канале доступен как `callvar.ann` — `tuple(sample_id, ann_vcf)`.

## Границы работы

**В работе:** `snp_sites`, `sample_snp_alleles`, `vcf_table`, `sample_snp_profiles`.

**Вне работы:** `vcf_storage`, per-sample `sample_software_provenance`, расхождения имён
файлов и колонок в `import_collection_batch.py` на стороне сайта, метаданные.

## Архитектура

Генерация разбита на два шага, потому что `run_batches.sh` гоняет 200k образцов
батчами по 500 с `--skip_final_reports`, а итоговые таблицы собираются один раз в конце.

```
main.nf (батч = 500 образцов)          batch_reports.nf (один раз в конце)
──────────────────────────────         ──────────────────────────────────
CALLVAR.out.ann ──> SNP_DB_SHARDS ──>  snp_db/shards/*.tsv.gz ──> SNP_DB_TABLES
                    (1 задача/батч)                                (склейка + дедуп)
                                                                        │
                                                          Reports/tb-platform/snp/
```

Одна задача Nextflow на батч, а не 200 000 задач на образец: внутри задачи параллелизм
даёт сам конвертер (`-j ${task.cpus}`).

### Компонент 1: `bin/vcf_to_snp_shards.py`

Основан на `/home/zerg/git/tb-platform/new_temp_data/vcf_to_snp_db_tsv_parallel_scan.py`.
Логика извлечения (`process_vcf_task`, `site_row_from_ann`, `ann_entries_for_alt`,
`gt_alt_indices`, `extract_metrics`, `parse_info`, `parse_format`) переносится **без
изменений** — она воспроизводит прод.

Что меняется:

- **удаляется** SQLite `site_map` и раздача `site_id`; на выходе `site_key`;
- **удаляется** запись `snp_sites.tsv`/`sample_snp_alleles.tsv`/`sample_snp_profiles.tsv`
  с числовыми ID, `schema_snp.sql`, `import_snp.sql`, `write_example_queries`;
- **удаляются** флаги `--snps-only`, `--pass-only`, `--min-dp`, `--min-af`,
  `--max-alt-alleles`, чтобы их нельзя было включить случайно и разойтись с продом;
- **добавляется** третий выход — строки `vcf_table`;
- **добавляется** дедупликация sites в пределах шарда с сохранением max(qual).

Интерфейс:

```
vcf_to_snp_shards.py --file-list vcfs.txt --tag <shard-tag> -o <outdir> -j <N>
```

Выход (без строки заголовка — заголовок добавляется только в итоговых файлах):

| Файл | Колонки |
|---|---|
| `sites.<tag>.tsv.gz` | `site_key, chrom, pos, ref, alt, qual, effect, impact, gene, gene_id, product_name, strand, feature_id, biotype, rank, hgvs_c, hgvs_p, cds_pos, aa_pos` |
| `alleles.<tag>.tsv.gz` | `sample_id, site_key, allele, dp, af, qual, filter` |
| `vcfrows.<tag>.tsv.gz` | `sample_id, pos, alt` |
| `samples.<tag>.txt` | список sample_id шарда (для манифеста) |

Дедупликация внутри шарда:

- `sites` — по `site_key`, оставляем строку с максимальным `qual`. Воркеры пишут
  per-file шарды во временный каталог, родительский процесс сливает их через словарь
  `site_key -> (row, qual)`. Оценка памяти: ~300k ключей × ~200 B ≈ 60 МБ.
- `alleles` — по `(sample_id, site_key)` в пределах файла (образец встречается в одном
  файле), шарды воркеров просто склеиваются.
- `vcfrows` — по `(sample_id, pos, alt)` в пределах файла.

Пустое значение в TSV = SQL NULL (`NULL ''` при COPY).

`sample_id` берётся из колонки образца в заголовке VCF, при её отсутствии — из имени
файла. Выбор источника зашит, отдельного флага нет. В tb-lite оба источника совпадают:
freebayes проставляет `SM:${meta.id}` из `@RG`.

### Компонент 2: `subworkflows/local/snp_db.nf` + `SNP_DB_SHARDS`

```groovy
process SNP_DB_SHARDS {
    label 'process_medium'
    container 'tb-lite/tb-platform-tables:1.1'   // python:3.12-slim, пересборка не нужна
    publishDir "${params.outdir}/snp_db/shards", mode: params.mode
    input:  path ann_vcfs
    output: path "sites.*.tsv.gz", path "alleles.*.tsv.gz",
            path "vcfrows.*.tsv.gz", path "samples.*.txt"
}
```

Вход — `CALLVAR.out.ann.map { _id, vcf -> vcf }.collect()`.

Тег шарда:

1. `params.batch_tag`, если задан — параметр уже есть в `nextflow.config`, и
   `run_batches.sh:413` уже передаёт `--batch_tag batch_${i}`, менять его не нужно;
2. иначе md5 отсортированного списка sample_id.

Оба варианта детерминированы: повторный прогон батча (`--resume-from`) перезаписывает
тот же файл, а не создаёт второй с дублями.

Подключение в `workflows/tblite.nf` — после `CALLVAR`, под гейтом `!params.skip_snp_db`.

Скрипт лежит в `bin/` и монтируется Nextflow'ом из `projectDir`, поэтому пересобирать
контейнер не требуется (как и для остальных `bin/*.py`).

### Компонент 3: `SNP_DB_TABLES` в `batch_reports.nf`

Собирает итог из шардов. Если каталог `snp_db/shards/` пуст — строит шарды сам из
`${params.outdir}/annotate_vcf/*/*.annotated.ann.vcf` через `.collate(500)`. Это
покрывает сценарий «пайплайн уже отработал, таблицы нужны задним числом» — например,
для готового прогона `/mnt/data1/tb-lite-runs/bc1286ed-5e2e-42db-a06c-0a5a02258824`.

Помимо вызова внутри `BATCH_REPORTS` есть отдельная точка входа
`batch_reports.nf -entry SNP_DB_ONLY` — только SNP-таблицы, без остальных отчётов.
Она нужна потому, что `BATCH_REPORTS` объявляет входные каналы с
`checkIfExists: true` и падает целиком, если в каталоге прогона пуст хотя бы один
из его глобов.

Склейка:

- **`snp_sites.tsv.gz`** — `sort -k1,1 -S 2G --parallel=${task.cpus} -T .` по всем
  `sites.*.tsv.gz`, затем awk оставляет на каждый `site_key` строку с максимальным
  `qual`. Это воспроизводит `ON CONFLICT … qual = max(qual)` из исходного скрипта.
- **`sample_snp_alleles.tsv.gz`** — конкатенация `alleles.*.tsv.gz`.
- **`vcf_table.tsv.gz`** — конкатенация `vcfrows.*.tsv.gz`.
- **`snp_db_manifest.tsv`** — число строк по каждому файлу, число образцов, число шардов.
- **`samples.txt`** — полный список sample_id.
- **`import_snp.sql`** — скрипт загрузки (ниже).

Заголовок дописывается только в итоговые файлы. Сжатие — `pigz -p ${task.cpus}`, если
доступен, иначе `gzip`; проверка через `command -v pigz`.

Публикация: `${params.outdir}/Reports/tb-platform/snp/`.

## Загрузка в БД: `import_snp.sql`

Генерируется `SNP_DB_TABLES` с подставленными абсолютными путями к трём `.tsv.gz`.
Запускается вручную через `psql`, **после** импорта `general.tsv` — на
`sample_snp_alleles.sample_id` висит FK на `general."ID"`.

```
1. TEMP staging stg_sites / stg_alleles / stg_vcf  <- \copy … FROM PROGRAM 'zcat …'
2. snp_sites          <- INSERT … ON CONFLICT (site_key) DO NOTHING
3. sample_snp_alleles <- INSERT … JOIN snp_sites USING (site_key)
                                  WHERE sample_id ∈ general
                                  ON CONFLICT (sample_id, site_id) DO NOTHING
4. vcf_table          <- DELETE USING (SELECT DISTINCT sample_id FROM stg_vcf),
                         затем INSERT (у таблицы нет уникального ключа)
5. sample_snp_profiles <- по 4 маскам, только для новых образцов, чанками по 5000,
                          ON CONFLICT (mask_key, sample_id) DO UPDATE
```

Профили — единственное, что нельзя посчитать в пайплайне: они хранят `site_id`, который
выдаёт sequence при вставке в `snp_sites`. Маски (`none`, `rlc_lowmap`, `dr_genes`,
`union`) берутся из `snp_mask_regions` в БД.

Канонизация site_id в профилях — как в `deploy/rebuild_snp_profiles_masked.sql`:

```sql
canon AS (SELECT site_id, min(site_id) OVER (PARTITION BY chrom,pos,ref,alt) AS cid
          FROM snp_sites
          WHERE char_length(ref)=1 AND char_length(alt)=1
            AND upper(ref) IN ('A','C','G','T','N') AND upper(alt) IN ('A','C','G','T','N'))
```

Существующие профили при этом не ломаются: sequence монотонна, новые `site_id` всегда
больше, поэтому канонический минимум для уже известных вариантов не меняется.

Агрегация по 200k образцов — это ~455M строк `sample_snp_alleles`, поэтому шаг 5
выполняется чанками по 5000 образцов (индекс `idx_sample_snp_alleles_sample`), а не
одним запросом.

Шаги 2–4 идут в одной транзакции, шаг 5 — отдельными транзакциями на чанк, чтобы
падение на профилях не откатывало загрузку данных.

## Идемпотентность

Повторный запуск импорта того же набора безопасен на каждом шаге:

| Таблица | Механизм |
|---|---|
| `snp_sites` | `ON CONFLICT (site_key) DO NOTHING` |
| `sample_snp_alleles` | `ON CONFLICT (sample_id, site_id) DO NOTHING` |
| `vcf_table` | DELETE по входящим `sample_id`, затем INSERT |
| `sample_snp_profiles` | `ON CONFLICT (mask_key, sample_id) DO UPDATE` |

Это же страхует от перекрытия шардов, если батч был перезапущен с изменившимся составом
образцов и появились два шарда с общими sample_id.

## Оценка объёмов (200 000 образцов)

Исходя из прод-соотношений (117k образцов → 130M строк `vcf_table`, 266M строк
`sample_snp_alleles`, 4.38M `snp_sites`):

| Артефакт | Строк | Raw | gz |
|---|---|---|---|
| `sample_snp_alleles.tsv.gz` | ~455M | ~32 ГБ | ~7 ГБ |
| `vcf_table.tsv.gz` | ~240M | ~5 ГБ | ~1 ГБ |
| `snp_sites.tsv.gz` | +2–4M к 4.38M | ~1 ГБ | ~0.2 ГБ |
| шарды `sites.*` до глобальной дедупликации | ~60M | ~9 ГБ | ~2 ГБ |

`sort` на шаге склейки `snp_sites` требует ~2× входа во временном каталоге задачи
(~18 ГБ) — задаётся через `-T .` и `-S 2G`.

## Тестирование

**Юнит, golden-file.** Три VCF из прогона
`/mnt/data1/tb-lite-runs/bc1286ed-5e2e-42db-a06c-0a5a02258824/results/annotate_vcf/`
(ERR4797591, ERR4797753, ERR4817364) прогоняются через `vcf_to_snp_shards.py`. Проверки:

- `site_key` для pos 1977 и 4013 равны прод-значениям из таблицы выше;
- число строк `vcfrows` для ERR4797591 = 1589 (все записи VCF, включая indel/complex/mnp);
- `product_name` и `strand` пустые во всех строках;
- `filter` пустой, `af` = 1;
- дедупликация: `site_key` в `sites.*` уникальны, `(sample_id, site_key)` в `alleles.*`
  уникальны, `(sample_id, pos, alt)` в `vcfrows.*` уникальны;
- у sites с одинаковым `site_key` от разных образцов остался максимальный `qual`.

**Интеграционный.** `import_snp.sql` прогоняется против **одноразового** контейнера
`postgres:17` со случайным именем: тест сам создаёт минимальную схему, заполняет её и
удаляет контейнер в teardown. К существующим базам тест не подключается — на
дев-машине в `tb_platform_db` лежит полная копия боевых данных, и незаскоупленный
`DELETE` там необратим. Тест выключен по умолчанию и включается переменной
`TB_LITE_SNP_DB_TEST=1`. Проверки:

- счётчики строк в `snp_sites` / `sample_snp_alleles` / `vcf_table` совпадают
  с `snp_db_manifest.tsv`;
- профили созданы для всех 4 масок, `snp_count` для маски `none` >= `snp_count`
  для маски `union`;
- повторный запуск того же `import_snp.sql` даёт 0 новых строк во всех четырёх таблицах.

**Проверка на живом пайплайне.** Только настоящие точки входа — у отдельного
тестового workflow в подкаталоге `projectDir` указывает не на корень репозитория,
и `${projectDir}/bin/...` не резолвится:

- `nextflow run batch_reports.nf -entry SNP_DB_ONLY --outdir <копия прогона
  bc1286ed-… без snp_db/>` — ветка автосборки шардов из `annotate_vcf/`;
- `nextflow run main.nf --input samplesheet.host.csv --batch_tag batch_1` на трёх
  образцах — связка `CALLVAR → SNP_DB`, и сверка `site_key` для pos 1977/4013
  с продовыми значениями на живом выводе пайплайна.

## Изменения в конфигурации

| Файл | Изменение |
|---|---|
| `nextflow.config` | `skip_snp_db = false` |
| `nextflow_schema.json` | описание `skip_snp_db` |
| `conf/modules.config` | `publishDir` для `SNP_DB_SHARDS` и `SNP_DB_TABLES` |
| `README.md` | раздел про SNP-таблицы и ручной импорт |
| `CHANGELOG.md` | запись о новых артефактах |

`run_batches.sh` не меняется: `params.batch_tag` уже передаётся.

Пересборка контейнеров не требуется: `SNP_DB_SHARDS` и `SNP_DB_TABLES` работают в уже
собранном `tb-lite/tb-platform-tables:1.1`, скрипты приезжают из `bin/`.

## Известные ограничения

- `sample_snp_profiles` строятся только на стороне БД — в пайплайне для этого нет
  `site_id`.
- Каталог `annotate_vcf/` при 200k образцов занимает ~282 ГБ в текущем виде (несжатый
  `.ann.vcf` 879 КБ + `genes.txt` 270 КБ + `snpEff_summary.html` 328 КБ на образец).
  Для SNP-таблиц публикация не нужна — `SNP_DB_SHARDS` читает файлы из work-каталога.
  Отключение публикации не входит в эту работу.
