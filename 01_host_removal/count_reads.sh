#!/usr/bin/env bash

FASTQ_DIR="/path/to/fastq_files"
OUTPUT="./read_counts.csv"

echo "sample,r1_reads,r2_reads,total_reads" > "$OUTPUT"

for R1 in "$FASTQ_DIR"/*_hostremoved_1.fq.gz; do
    SAMPLE=$(basename "$R1" _hostremoved_1.fq.gz)
    R2="${R1/_hostremoved_1.fq.gz/_hostremoved_2.fq.gz}"
    R1_READS=$(zcat "$R1" | awk 'END {print NR/4}')
    R2_READS=$(zcat "$R2" | awk 'END {print NR/4}')
    echo "${SAMPLE},${R1_READS},${R2_READS},$(( R1_READS + R2_READS ))" >> "$OUTPUT"
done