#!/usr/bin/env bash
# Usage: count_reads.sh <output_tsv> <accessions_file> <clean_bam1> [clean_bam2 ...]
# For each sample's final, fully-cleaned BAM (remove_flagged_regions output),
# and for each accession/contig listed in <accessions_file>, counts unique
# read names with duplicates included, and with duplicates excluded
# (SAM flag 1024).
set -euo pipefail
OUT="$1"; ACCESSIONS="$2"; shift 2

printf "sample\taccession\treads_with_dups\treads_without_dups\n" > "$OUT"

for BAM in "$@"; do
    SAMPLE=$(basename "$BAM" .bam)

    # index if needed
    [[ -f "${BAM}.bai" || -f "${BAM%.bam}.bai" ]] || samtools index "$BAM"

    for ACC in $(cat "$ACCESSIONS"); do
        WITH=$(samtools view "$BAM" "$ACC" 2>/dev/null | cut -f1 | sort -u | wc -l)
        NODUP=$(samtools view -F 1024 "$BAM" "$ACC" 2>/dev/null | cut -f1 | sort -u | wc -l)
        printf "%s\t%s\t%s\t%s\n" "$SAMPLE" "$ACC" "$WITH" "$NODUP" >> "$OUT"
    done
done
