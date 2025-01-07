### Inputs 
#### data/

- SNP_table_noresis.txt (SNP table WITHOUT resistance mutations)

- original alignment (check that ref and lineage0 are included)

- original tree (converted to relative distances and with cleaned name)



- mysplits_column.txt


### Pipeline
#### homoplasy_finder/

1. Sort alignment

```bash
seqkit sort -N alignment
```

2. Split the alignment into 200 parts

```bash
splitalignment.py --aln alignment.fas --fnr 200 --dir split_results
```

3. Launch `treetime` for each alignment chunk

```bash
parallel -j 3 treetime ancestral --aln {} --tree ../data/tree_RELATIVE_CLEANED.nwk --outdir ancestral_results/{.} ::: split_results/*.fas
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