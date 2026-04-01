#!/usr/bin/env bash
set -euo pipefail

# add_unclassified.sh
# Adds "Unclassified" row to Bracken output using Kraken2 report
# Usage: add_unclassified.sh -o output.tsv bracken.tsv kraken_report.txt

OUTPUT=""
while getopts "o:" opt; do
    case $opt in
        o) OUTPUT="$OPTARG" ;;
        *) echo "Usage: $0 -o output.tsv bracken.tsv kraken_report" >&2; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

BRACKEN_TSV="$1"
KRAKEN_REPORT="$2"

# Get unclassified count from kraken report (first field of the line starting with unclassified)
UNCLASSIFIED=$(awk '$4 == "U" { print $2; exit }' "$KRAKEN_REPORT")
if [ -z "$UNCLASSIFIED" ]; then
    UNCLASSIFIED=0
fi

# Get total reads from kraken report
TOTAL=$(awk 'NR==1 { gsub(/ /, "", $2); print $2 }' "$KRAKEN_REPORT")
if [ -z "$TOTAL" ] || [ "$TOTAL" -eq 0 ]; then
    TOTAL=1
fi

# Calculate fraction
FRAC=$(awk "BEGIN { printf \"%.6f\", ${UNCLASSIFIED}/${TOTAL} }")

# Copy bracken output and append unclassified row
cp "$BRACKEN_TSV" "$OUTPUT"
echo -e "Unclassified\t0\tS\t${UNCLASSIFIED}\t0\t${UNCLASSIFIED}\t${FRAC}" >> "$OUTPUT"
