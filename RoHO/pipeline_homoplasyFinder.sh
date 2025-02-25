#!/bin/bash

# 1. Sort alignment
seqkit sort -N ../data/alignment.fas

# 2. Split the alignment into 200 parts
splitalignment.py --aln ../data/alignment.fas --fnr 200 --dir ../data/split_results

# 3. Launch `treetime` for each alignment chunk
parallel -j 10 treetime ancestral --aln {} --tree ../data/tree.nwk --outdir ../data/ancestral_results/{/.} ::: ../data/split_results/*.fas

# 4. Merge and parse mutation positions in partial alignments
Rscript mergeAncestralResults.R

# 5. Select SNP table mutations
Rscript filterSNPtable.R

# 6. Filter SNPs that appear more than once in the tree
Rscript filterHomoplasySNPtable.R

# 7. Find homoplasies
Rscript homoplasyFinder.R

# 8. Add gene names to found homoplasy
Rscript addGeneToHomoplasy.R