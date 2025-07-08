#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# build-containers.sh  ─  собрать или скачать все Singularity‑образы проекта
# -----------------------------------------------------------------------------
#
# • Строит *.sif из *.def (Singularity recipes) в каталоге containers/.
# • Скачивает *.sif, перечисленные в containers/pull-images.txt.
#   Формат файла: <sif‑name> <URI>, например:
#       fastqc_latest.sif  docker://staphb/fastqc:latest
#       multiqc.sif        oras://ghcr.io/org/multiqc:1.15
# • Пропускает сборку/скачивание, если .sif свежее рецепта или уже существует
#   (перестраивать/перекачивать можно опцией --force).
# • Работаeт как с Singularity, так и c Apptainer.
# -----------------------------------------------------------------------------
set -euo pipefail

###############################################################################
# ПАРАМЕТРЫ КОМАНДНОЙ СТРОКИ                                                  #
###############################################################################
FORCE=0         # 1 → пересобрать / перекачать даже если образ актуален
FAKEROOT=0      # 1 → использовать --fakeroot при сборке .def

usage() {
  cat <<EOF
Usage: ${0##*/} [options]
Options:
  -f, --force        rebuild / repull even if image exists and is up‑to‑date
  -r, --fakeroot     use --fakeroot flag when building from .def
  -h, --help         show this help and exit
EOF
}

while [[ $# -gt 0 ]]; do
  case $1 in
    -f|--force) FORCE=1 ; shift ;;
    -r|--fakeroot) FAKEROOT=1 ; shift ;;
    -h|--help) usage ; exit 0 ;;
    *) echo "Unknown option: $1" >&2 ; usage ; exit 1 ;;
  esac
done

###############################################################################
# ОПРЕДЕЛЯЕМ БИНАРНИК SINGULARITY / APPTAINER                                 #
###############################################################################
get_singularity() {
  if command -v singularity &>/dev/null; then echo singularity; return 0; fi
  if command -v apptainer   &>/dev/null; then echo apptainer;   return 0; fi
  echo "Error: neither singularity nor apptainer found in \$PATH" >&2; exit 1
}
SINGULARITY=$(get_singularity)

###############################################################################
# КАТАЛОГИ И ФАЙЛЫ                                                            #
###############################################################################
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CONTAINERS_DIR="$SCRIPT_DIR/containers/def"
PULL_MANIFEST="$CONTAINERS_DIR/pull-images.txt"

if [[ ! -d $CONTAINERS_DIR ]]; then
  echo "Error: $CONTAINERS_DIR not found" >&2
  exit 1
fi

###############################################################################
# ФУНКЦИЯ: СБОРКА .def → .sif                                                #
###############################################################################
build_defs() {
  shopt -s nullglob
  local defs=("$CONTAINERS_DIR"/*.def)
  shopt -u nullglob

  for def in "${defs[@]}"; do
    local base=$(basename "$def" .def)
    local sif="$CONTAINERS_DIR/$base.sif"

    local need_build=0
    if [[ $FORCE -eq 1 ]]; then
      need_build=1
    elif [[ ! -f $sif ]]; then
      need_build=1
    elif [[ $def -nt $sif ]]; then
      need_build=1
    fi

    if [[ $need_build -eq 1 ]]; then
      echo "[+] Building $sif ← $def"
      local cmd=($SINGULARITY build)
      [[ $FAKEROOT -eq 1 ]] && cmd+=(--fakeroot)
      cmd+=("$sif" "$def")
      "${cmd[@]}"
    else
      echo "[=] Skipping $sif (up–to‑date)"
    fi
  done
}

###############################################################################
# ФУНКЦИЯ: ЗАГРУЗКА .sif согласно pull-images.txt                            #
###############################################################################
pull_images() {
  [[ -f $PULL_MANIFEST ]] || return 0   # ничего качать

  while read -r line; do
    # пропускаем пустые строки и комментарии
    [[ -z "$line" || "${line:0:1}" == "#" ]] && continue

    # парсим: первое поле — имя файла, всё остальное — URI
    local sif uri
    sif=$(echo "$line" | awk '{print $1}')
    uri=$(echo "$line" | cut -d' ' -f2-)

    [[ -z $sif || -z $uri ]] && {
      echo "Warning: malformed line in $PULL_MANIFEST: $line" >&2; continue; }

    local dest="$CONTAINERS_DIR/$sif"
    local need_pull=0

    if [[ $FORCE -eq 1 ]]; then
      need_pull=1
    elif [[ ! -f $dest ]]; then
      need_pull=1
    fi

    if [[ $need_pull -eq 1 ]]; then
      echo "[+] Pulling $dest ← $uri"
      "$SINGULARITY" pull --name "$dest" "$uri"
    else
      echo "[=] Skipping $dest (already exists)"
    fi

  done < "$PULL_MANIFEST"
}

###############################################################################
# ЗАПУСК                                                                       #
###############################################################################
build_defs
pull_images

echo "✓ All container images are ready."
