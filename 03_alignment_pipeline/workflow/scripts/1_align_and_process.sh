#!/usr/bin/env bash
# Usage: align_and_process.sh <fq1> <fq2> <ref> <threads> <sample> <out_bam>
set -euo pipefail
FQ1="$1"; FQ2="$2"; REF="$3"; THREADS="$4"; SAMPLE="$5"; OUT_BAM="$6"

OUTDIR=$(dirname "$OUT_BAM")
SAM="$OUTDIR/${SAMPLE}.sam"
HEADER="$OUTDIR/${SAMPLE}_header.sam"
READS="$OUTDIR/${SAMPLE}_reads.sam"
NAMESORT="$OUTDIR/${SAMPLE}.namesort.bam"
FIXMATE="$OUTDIR/${SAMPLE}.fixmate.bam"
SORTED="$OUTDIR/${SAMPLE}.sorted.bam"

bwa-mem2 mem -t "$THREADS" "$REF" "$FQ1" "$FQ2" > "$SAM"

samtools view "$SAM" \
    | awk -v s="$SAMPLE" 'BEGIN{OFS="\t"} /^@/{print;next} {$1=s"_"$1; print}' > "$READS"
samtools view -H "$SAM" > "$HEADER"

cat "$HEADER" "$READS" \
    | samtools sort -n -@ "$THREADS" -o "$NAMESORT"

samtools fixmate -m -@ "$THREADS" "$NAMESORT" "$FIXMATE"
samtools sort -@ "$THREADS" -o "$SORTED" "$FIXMATE"
samtools markdup -@ "$THREADS" -s "$SORTED" "$OUT_BAM"
samtools index "$OUT_BAM"

rm -f "$SAM" "$HEADER" "$READS" "$NAMESORT" "$FIXMATE" "$SORTED"
