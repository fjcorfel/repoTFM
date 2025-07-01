library(treeio)
library(ape)

# Load tree
tree <- read.tree("../data/tree.nwk")


# Convert absolute distances to relative distances
# Alignment lenght = 4798 (extracted from mysplits.txt with 200 fragments)
tree$edge.length <- tree$edge.length / 4798


# Write new nwk tree with relative distances
write.tree(tree, file = "tree_RELATIVE.nwk")
