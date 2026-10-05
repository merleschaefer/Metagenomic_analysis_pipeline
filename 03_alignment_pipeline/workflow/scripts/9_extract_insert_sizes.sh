#!/usr/bin/env bash
# Usage: extract_insert_sizes.sh <output_txt> <clean_bam1> [clean_bam2 ...]
set -euo pipefail
OUT="$1"; shift

echo -e "size\tsample" > "$OUT"
for BAM in "$@"; do
    SAMPLE=$(basename "$BAM" .clean.bam)
    samtools view -f 0x2 -F 1024 "$BAM" | awk -v s="$SAMPLE" '$9 > 0 {print $9"\t"s}' >> "$OUT"
done
