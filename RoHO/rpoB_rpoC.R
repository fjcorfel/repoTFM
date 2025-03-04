library(dplyr)
library(phangorn)
library(phytools)

load("../../data_RoHO/homoplasy_nodes_annotated.rda")
tree <- ape::as.phylo(treeio::read.beast("../../data_RoHO/annotated_tree_resis.nexus")) 

rpoB_nodes <- homoplasy_nodes %>%
  filter(synonym == "rpoB") %>%
  mutate(has_rpoC = FALSE,
         rpoC_node = vector("list", n()))

rpoC_nodes <- homoplasy_nodes %>%
  filter(synonym == "rpoC")


# Check rpoB nodes that have rpoC (ancestors and descendants)
for (n_node in seq_along(rpoB_nodes$node)) {
  node = rpoB_nodes$node[n_node]
  mutation = rpoB_nodes$mutation[n_node]
  
  # Check for rpoC in node descendants
  node_descendants <- unlist(phangorn::Descendants(tree, node, "all")) 
  
  matching_nodes <- node_descendants[node_descendants %in% rpoC_nodes$node]
  

  if (length(matching_nodes) > 0) {
    rpoB_nodes$has_rpoC[n_node] <- TRUE
    rpoB_nodes$rpoC_node[[n_node]] <- matching_nodes
  }
}

rpoB_nodes_rpoC <- rpoB_nodes %>%
  filter(has_rpoC)

rpoC_nodes <- unlist(rpoB_nodes_rpoC$rpoC_node)
rpoC_nodes_with_rpoB <- data.frame(rpoC_node = rpoC_nodes)

# Get mut alleles from a given node (check node tips )
get_mut_alleles <- function(node) {
  tips <- unlist(phangorn::Descendants(tree, node, "tips"))
  
  return(length(tips))
}

# Get wt alleles from a given node (check sister tips)
get_wt_alleles <- function(node) {
  sister <- phangorn::Siblings(tree, node)
  sister_tips <- unlist(phangorn::Descendants(tree, sister, "tips"))
  
  return(length(sister_tips))
}

rpoC_nodes_with_rpoB <- rpoC_nodes_with_rpoB %>%
  rowwise() %>%
  mutate(mutation = homoplasy_nodes$mutation[match(rpoC_node, homoplasy_nodes$node)]) %>%
  mutate(
    n_mut_alleles = get_mut_alleles(rpoC_node)
  ) %>%
  mutate(
    n_wt_alleles = get_wt_alleles(rpoC_node)
  ) %>%
  mutate(
    RoHO = n_mut_alleles / n_wt_alleles,
    sister_RoHO = n_wt_alleles / n_mut_alleles
  )


t.test(rpoC_nodes_with_rpoB$RoHO, rpoC_nodes_with_rpoB$sister_RoHO, alternative = "greater")$p.value
wilcox.test(rpoC_nodes_with_rpoB$RoHO, rpoC_nodes_with_rpoB$sister_RoHO, alternative = "greater")$p.value

rpoB_nodes_without_rpoC <- rpoB_nodes %>%
  filter(!has_rpoC)

rpoB_nodes_without_rpoC <- rpoB_nodes_without_rpoC %>%
  rowwise() %>%
  mutate(
    n_mut_alleles = get_mut_alleles(node)
  ) %>%
  mutate(
    n_wt_alleles = get_wt_alleles(node)
  ) %>%
  mutate(
    RoHO = n_mut_alleles / n_wt_alleles,
    sister_RoHO = n_wt_alleles / n_mut_alleles
  )

t.test(rpoB_nodes_without_rpoC$RoHO, rpoB_nodes_without_rpoC$sister_RoHO, alternative = "less")$p.value
wilcox.test(rpoB_nodes_without_rpoC$RoHO, rpoB_nodes_without_rpoC$sister_RoHO, alternative = "less")$p.value
