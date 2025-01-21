### Inputs 
#### data/

- SNP_table_noresis.txt (SNP table WITHOUT resistance mutations)

- alignment.fas (check that ref and lineage0 are included)

- tree_relative_cleaned.nwk (converted to relative distances and with cleaned name)

- mysplits_column.txt


### Pipeline
#### homoplasy_finder/

1. Sort alignment

```bash
seqkit sort -N ../data/alignment.fas
```

2. Split the alignment into 200 parts

```bash
splitalignment.py --aln ../data/alignment.fas --fnr 200 --dir ../data/split_results
```

3. Launch `treetime` for each alignment chunk

```bash
parallel -j 3 treetime ancestral --aln {} --tree ../data/tree_relative_cleaned.nwk --outdir ../data/ancestral_results/{.} ::: ../data/split_results/*.fas
```

4. Parse mutation positions in partial alignments to genome positions

```bash
Rscript mergeAncestralResults.R
```

#### WIP
5. Find homoplasies from mutation data 

```bash
Rscript homoplasyFinder.R
```