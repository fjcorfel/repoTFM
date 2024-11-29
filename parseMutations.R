library(treeio)
library(ape)
library(dplyr)

# Read SNP table -> df
snpTable <- read.table("./data/SNP_table_noresis.txt", sep = "\t", header = TRUE)

# Read splits -> vector
splits <- readLines("./data/mysplits.txt")
splits <- gsub("\\[|\\]", "", splits)
splits <- as.numeric(unlist(strsplit(splits, ", ")))



# For each element of the split, open a tree




# Read tree from file -> phylo object
#tree <- read.beast(file = "annotated_trees/annotated_tree_001.nexus")
# Tree tibble conversion
#tibble_tree <- as_tibble(tree)
