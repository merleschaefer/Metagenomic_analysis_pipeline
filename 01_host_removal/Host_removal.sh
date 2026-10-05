#!/bin/bash
set -euo pipefail

# ── Config ────────────────────────────────────────────────────────────────────
FQ_DIR="/path/to/FASTQ_FILES" 
HOST_REF="/path/to/host_reference_fasta"  
OUTDIR="./results" 
THREADS=16
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p "$OUTDIR/bam" "$OUTDIR/fastq"

[ -f "${HOST_REF}.0123" ] || bwa-mem2 index "$HOST_REF"

for FQ1 in "$FQ_DIR"/*_R1_001.fastq.gz; do # adjust to match endings of your FastQ files!
    SAMPLE=$(basename "$FQ1" _R1_001.fastq.gz) # adjust to match endings of your FastQ files!
    FQ2="$FQ_DIR/${SAMPLE}_R2_001.fastq.gz" # adjust to match endings of your FastQ files!

    BAM_NAMESORT="$OUTDIR/bam/${SAMPLE}.namesort.bam"
    OUT_FQ1="$OUTDIR/fastq/${SAMPLE}_hostremoved_1.fq.gz"
    OUT_FQ2="$OUTDIR/fastq/${SAMPLE}_hostremoved_2.fq.gz"

    echo "==> Host-removing $SAMPLE"

    # Align to host, pipe directly to name-sort
    bwa-mem2 mem -t "$THREADS" "$HOST_REF" "$FQ1" "$FQ2" \
        | samtools sort -n -@ "$THREADS" -o "$BAM_NAMESORT"

    # -f 12: include reads when both reads are unmapped (neither maps to host)
    # -F 256: exclude secondary alignments to host
    samtools fastq \
        -@ "$THREADS" \
        -f 12 -F 256 \
        -1 "$OUT_FQ1" \
        -2 "$OUT_FQ2" \
        "$BAM_NAMESORT"

    rm "$BAM_NAMESORT"

    echo "==> Done: $OUT_FQ1 / $OUT_FQ2"
done