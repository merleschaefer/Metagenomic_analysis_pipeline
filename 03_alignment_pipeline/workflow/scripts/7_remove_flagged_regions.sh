#!/usr/bin/env bash
# Usage: remove_flagged_regions.sh <filtered_bam> <flagged_regions.bed> <output_clean_bam>
set -euo pipefail
BAM="$1"; BED="$2"; OUT="$3"

BAD="${OUT%.bam}.bad"
GOOD="${OUT%.bam}.good"

bedtools intersect -abam "$BAM" -b "$BED" -u \
    | samtools view - | cut -f1 | sort -u > "$BAD"

samtools view "$BAM" | cut -f1 | sort -u \
    | comm -23 - "$BAD" > "$GOOD"

samtools view -N "$GOOD" -b "$BAM" > "$OUT"

rm -f "$BAD" "$GOOD"
