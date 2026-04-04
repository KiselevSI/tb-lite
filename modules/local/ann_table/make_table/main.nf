process MAKE_TABLE {
  tag        "SNPEFF_table"
  label 'process_low'
  publishDir "${params.outdir}/Reports/snp_matrix", mode: params.mode

  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'docker://tb-lite/ann-table:1.0' :
      'tb-lite/ann-table:1.0' }"

  input:
    tuple path(vcf), path(tbi)

  output:
    path "ANNOTATION_TABLE.tsv"

  script:
  """
  SNPSIFT_BIN=\$(command -v SnpSift || command -v snpSift || true)
  [[ -n "\$SNPSIFT_BIN" ]] || { echo "Neither SnpSift nor snpSift found in PATH" >&2; exit 127; }

  # 1) Левый блок c шапкой от SnpSift (добавлены *_LEN)
  "\$SNPSIFT_BIN" extractFields -e "." -s "," "$vcf" \
    CHROM POS REF ALT QUAL \
    "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" \
    "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" \
    "ANN[*].HGVS_C" "ANN[*].HGVS_P" \
    "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" "ANN[*].AA_LEN" \
    > ann_core.withhdr.tsv

  # 2) Склеиваем POS/LEN и выбрасываем *_LEN (без функции, совместимо с busybox/mawk)
awk -F'\t' -v OFS='\t' '
NR==1{
  # индексы колонок по имени
  for(i=1;i<=NF;i++){ idx[\$i]=i }
  cp=idx["ANN[*].CDS_POS"]; cl=idx["ANN[*].CDS_LEN"];
  ap=idx["ANN[*].AA_POS"];  al=idx["ANN[*].AA_LEN"];

  # какие колонки оставить (исключаем *_LEN)
  keep_cnt=0
  for(i=1;i<=NF;i++){
    if(i==cl || i==al) continue
    keep[++keep_cnt]=i
  }

  # печатаем заголовок (без *_LEN)
  out=\$(keep[1])
  for(k=2;k<=keep_cnt;k++){ out = out OFS \$(keep[k]) }
  print out
  next
}
{
  # helper: склейка списков pos,len -> "pos/len"
  if (cp && cl){
    n1=split(\$(cp), pa, ","); n2=split(\$(cl), pl, ",")
    n=(n1>n2?n1:n2); res=""
    for(i=1;i<=n;i++){
      pv=pa[i]; lv=pl[i]
      if(pv=="" || pv=="." || pv=="-1" || lv=="" || lv=="."){ v="." } else { v=pv "/" lv }
      res = (res=="" ? v : res "," v)
    }
    \$(cp)=res
  }
  if (ap && al){
    n1=split(\$(ap), pa, ","); n2=split(\$(al), pl, ",")
    n=(n1>n2?n1:n2); res=""
    for(i=1;i<=n;i++){
      pv=pa[i]; lv=pl[i]
      if(pv=="" || pv=="." || pv=="-1" || lv=="" || lv=="." ){ v="." } else { v=pv "/" lv }
      res = (res=="" ? v : res "," v)
    }
    \$(ap)=res
  }

  # печатаем только выбранные колонки
  out=\$(keep[1])
  for(k=2;k<=keep_cnt;k++){ out = out OFS \$(keep[k]) }
  print out
}
' ann_core.withhdr.tsv > ann_core.fixed.withhdr.tsv

  # 3) Общий заголовок: левая шапка + имена сэмплов
  LEFT_HDR=\$(head -n1 ann_core.fixed.withhdr.tsv)
  SAMPLES_TSV=\$(bcftools query -l "$vcf" | paste -sd \$'\\t' -)
  printf "%s\\t%s\\n" "\$LEFT_HDR" "\$SAMPLES_TSV" > header.tsv

  # 4) Данные слева без шапки
  tail -n +2 ann_core.fixed.withhdr.tsv > ann_core.tsv

  # 5) Правый блок: только генотипы буквами, без шапки
  SAMPLES_CSV=\$(bcftools query -l "$vcf" | paste -sd, -)
  bcftools query -s "\$SAMPLES_CSV" -f '[\\t%TGT]\\n' "$vcf" | sed 's/^\\t//' > genotypes_only.tsv

  # 6) Склейка в финальную таблицу
  cat header.tsv > ANNOTATION_TABLE.tsv
  paste ann_core.tsv genotypes_only.tsv >> ANNOTATION_TABLE.tsv
  """
}
