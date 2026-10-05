#!/bin/bash

# ── Paths ────────────────────────────────────────────────
INPUT="./samplesheet.csv"
DATABASES="./databasesheet.csv"
OUTDIR="./results"

# ── Run pipeline ─────────────────────────────────────────
nextflow run nf-core/taxprofiler \
    -r 2.0.1 \
    --input "$INPUT" \
    --databases "$DATABASES" \
    --outdir "$OUTDIR" \
    -profile docker \
    \
    `# ── Pre-processing ──────────────────────────────────` \
    --perform_shortread_qc \
    --shortread_qc_minlength 130 \
    --shortread_qc_dedup \
    \
    `# ── Classification ──────────────────────────────────` \
    --run_kraken2 \
    --run_bracken \
    --bracken_save_intermediatekraken2 \
    --run_profile_standardisation \
    \
    -resume