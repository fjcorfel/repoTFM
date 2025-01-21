library(treeio)
library(ape)

# Load tree (phylo objetc) and convert it to tibble
tree <- treeio::read.beast.newick("tree.nwk")
tibble_tree <- as_tibble(tree)

# Convert absolute distances to relative distances
# Alignment lenght = 4798 (extracted from mysplits.txt with 200 fragments)
tibble_tree$branch.length <- tibble_tree$branch.length / 4798

# Revert tibble to phylo object
relative_tree <- treeio::as.phylo(tibble_tree)

# Write new nwk tree with relative distances
treeio::write.tree(relative_tree, file = "tree_RELATIVE.nwk")
