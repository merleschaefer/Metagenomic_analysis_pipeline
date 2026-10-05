#!/bin/bash

# Directory containing your fastq files
FASTQ_DIR="/path/to/hostremoved_FastQ_files"
# Output samplesheet
OUTPUT="./samplesheet.csv"

# Write header
echo "sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta" > "$OUTPUT"

# Loop over R1 files (one per sample)
for R1 in "$FASTQ_DIR"/*_hostremoved_1.fq.gz; do
    # Get sample name
    BASENAME=$(basename "$R1" _hostremoved_1.fq.gz)
    
    # Build R2 path
    R2="$FASTQ_DIR/${BASENAME}_hostremoved_2.fq.gz"
    
    # Write row — adjust ILLUMINA and run accession as needed
    echo "${BASENAME},${BASENAME},ILLUMINA,${R1},${R2}," >> "$OUTPUT"
done

echo "Samplesheet written to $OUTPUT"