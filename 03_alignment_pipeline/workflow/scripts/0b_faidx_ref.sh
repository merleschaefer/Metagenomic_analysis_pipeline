#!/usr/bin/env bash
# Usage: faidx_ref.sh <reference.fasta>
set -euo pipefail
REF="$1"
samtools faidx "$REF"
