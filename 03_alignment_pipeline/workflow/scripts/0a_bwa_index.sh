#!/usr/bin/env bash
# Usage: bwa_index.sh <reference.fasta>
set -euo pipefail
REF="$1"
bwa-mem2 index "$REF"
