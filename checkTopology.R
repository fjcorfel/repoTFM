library(treeio)
library(data.table)
library(phytools)

original_tree <- treeio::read.newick("./global_rescatado_RELATIVE_CLEANED.nwk")
original_tibble_tree <- as_tibble(original_tree)

annotated_tree <- read.beast("./annotated_tree.nexus")
annotated_tree <- as.phylo(annotated_tree)
annotated_tibble_tree <- as_tibble(annotated_tree)

# node_55050
original_node <- 39685
length(phytools::getDescendants(original_tree, original_node))
annotated_node <- 42855
length(phytools::getDescendants(annotated_tree, annotated_node))


# node_40053
original_node <- 61614
length(phytools::getDescendants(original_tree, original_node))
annotated_node <- 49399
length(phytools::getDescendants(annotated_tree, annotated_node))


# node_56063
original_node <- 42972
original_descendants <- phytools::getDescendants(original_tree, original_node)
annotated_node <- 40553
annotated_descendants <- phytools::getDescendants(annotated_tree, annotated_node)

x <- original_tibble_tree[which(original_tibble_tree$node %in% original_descendants), ]
y <- annotated_tibble_tree[which(annotated_tibble_tree$node %in% annotated_descendants), ]

x <- x$label
y <- y$label

setdiff(x, y)
setdiff(y, x)
