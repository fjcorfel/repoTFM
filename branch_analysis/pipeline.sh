#!/bin/bash

BASE_DIR="../data/POPULATION_BASED"

DATASETS=($(find "$BASE_DIR" -mindepth 1 -maxdepth 1 -type d -exec basename {} \;))

echo "List of datasets: ${DATASETS[@]}"
echo -e "\n"

echo "1. Running TreeTime ancestral reconstruction for all datasets..."
./run_treetime.sh || { echo "Error in TreeTime ancestral reconstruction"; exit 1; }
Rscript nexus_to_csv.R || { echo "Error in nexus_to_csv.R"; exit 1; }
echo -e "\n\nTreeTime ancestral reconstruction has been completed.\n\n"

# Process each dataset
for dataset in "${DATASETS[@]}"; do
    echo -e "\n"
    echo "----------------------------------------"
    echo "Processing dataset: $dataset"
    echo "----------------------------------------"
    echo -e "\n"

    echo "2. Cleaning $dataset TreeTime results..."
    python clean_annotated_tree.py --dataset "$dataset" || { echo "Error in clean_annotated_tree.py - $dataset"; continue; }

    echo "3. Counting SNPs in $dataset dataset..."
    python count_SNPs.py --dataset "$dataset" || { echo "Error in count_SNPs.py - $dataset"; continue; }

    echo "4. Performing branch length analysis for $dataset dataset..."
    Rscript branch_length_analysis.R "$dataset" || { echo "Error in branch_length_analysis.R - $dataset"; continue; }

    echo -e "Processing of $dataset dataset completed.\n\n"
done
