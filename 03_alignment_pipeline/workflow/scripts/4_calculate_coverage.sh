#!/usr/bin/env bash
# Usage: calculate_coverage.sh <merged_filtered_bam> <accessions_file> <output_tsv>
set -euo pipefail
BAM="$1"; ACCESSIONS="$2"; OUT="$3"

samtools view -b -F 1024 "$BAM" $(cat "$ACCESSIONS" | tr '\n' ' ') \
    | samtools depth -a - > "$OUT"
