library(treeio)
library(ape)
library(dplyr)

# Read mysplits.txt
splits <- readLines("mysplits.txt")

# Read tree from file -> phylo object
tree <- read.beast(file = "annotated_trees/annotated_tree_001.nexus")

# Tree tibble conversion
tibble_tree <- as_tibble(tree)
