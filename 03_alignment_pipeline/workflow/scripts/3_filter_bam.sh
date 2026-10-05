#!/usr/bin/env bash
# Usage: filter_bam.sh <input_bam> <output_bam> <filter_expr> <min_mapq> <require_flag>
set -euo pipefail
IN="$1"; OUT="$2"; EXPR="$3"; MAPQ="$4"; FLAG="$5"

samtools view \
    -e "$EXPR" \
    -q "$MAPQ" \
    -f "$FLAG" \
    -b "$IN" \
    -o "$OUT"
samtools index "$OUT"
