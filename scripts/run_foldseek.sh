#!/bin/bash
# Run Foldseek all-vs-all structural search on all PDB files.
# Uses the official Steinegger lab Docker image (no build required).
# Output: ~/pen-stack/data/foldseek/results.m8
#
# Usage: bash scripts/run_foldseek.sh
# Runtime: ~30-60 min on 24-core VM for 2,284 PDB structures
#
# Column order in results.m8 (set by --format-output below):
#   query target fident alnlen mismatch gapopen qstart qend tstart tend evalue bits tmscore
#   (tmscore is the TM-score from TMalign, range 0-1, >0.5 = same fold)

set -euo pipefail

DATA=/home/anees_22phd0670/pen-stack/data
PDB_DIR=$DATA/raw/pdb
OUT_DIR=$DATA/foldseek
TMP_DIR=/tmp/foldseek_tmp

mkdir -p "$OUT_DIR" "$TMP_DIR"

echo "[$(date)] Starting Foldseek all-vs-all search on $PDB_DIR"
echo "[$(date)] Output -> $OUT_DIR/results.m8"

docker run --rm \
    -v "$DATA:/data" \
    -v "$TMP_DIR:/tmp/foldseek_tmp" \
    ghcr.io/steineggerlab/foldseek:latest \
    easy-search \
        /data/raw/pdb \
        /data/raw/pdb \
        /data/foldseek/results.m8 \
        /tmp/foldseek_tmp \
        --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tmscore" \
        --tmscore-threshold 0.0 \
        --max-seqs 50 \
        --threads 20 \
        -v 2

echo "[$(date)] Foldseek done."
echo "Lines in results: $(wc -l < $OUT_DIR/results.m8)"
echo "Next step: python3 scripts/add_structural_edges.py"
