#!/usr/bin/env bash
# Usage: merge_bam.sh <threads> <output_bam> <input_bam1> [input_bam2 ...]
set -euo pipefail
THREADS="$1"; OUT="$2"; shift 2

samtools merge -f -@ "$THREADS" "$OUT" "$@"
samtools index "$OUT"
