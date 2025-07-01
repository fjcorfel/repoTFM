#!/bin/bash

N_CORES=8

# Main directory with population based datasets
BASE_DIR="../data/POPULATION_BASED"

# Find alignment file (.fas) inside directories
# e.g. "POPULATION_BASED/Valencia/*.fas"
find "$BASE_DIR" -mindepth 2 -type f -name "*.fas" | parallel -j "$N_CORES" '
    ALIGNMENT={}
    DIRNAME=$(dirname "$ALIGNMENT")

    OUT_DIR="$DIRNAME/treetime_results/"
    mkdir -p "$OUT_DIR"

    treetime ancestral \
        --aln "$ALIGNMENT" \
        --tree "$DIRNAME/tree.nwk" \
        --outdir "$OUT_DIR"
    '

echo -e "\n\nTreeTime ancestral reconstruction has been completed.\n\n"
