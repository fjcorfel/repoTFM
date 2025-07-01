library(tidyverse)
library(phangorn)
library(data.table)

# Load tree (nwk)
DATASET <- "vietnam"
tree <- ape::as.phylo(treeio::read.beast(paste0("../data/", DATASET, "/", "annotated_tree.nexus"))) 
annotated_tree <- fread(paste0("../data/", DATASET, "/", "annotated_tree_cleaned.csv")) %>%
  as_tibble()


get_ancestors <- function(node) {
  ancestors <- phangorn::Ancestors(tree, node)
  return(ancestors)
}  

get_descendants <- function(node) {
  descendants <- phangorn::Descendants(tree, node, "all")
  return(descendants)
}  

get_n_synonym_mutations <- function(node) {
  mutation_idx <- which(annotated_tree$node == node)
  if (length(mutation_idx) > 0) {
    return(annotated_tree$n_synonym_mutations[mutation_idx])
  } else return(0)
}

accumulated_mutations_root <- sapply(annotated_tree$node, function(n) {
  node_mutations <- get_n_synonym_mutations(n)
  ancestors <- get_ancestors(n)
  ancestors_mutations <- sum(sapply(ancestors, function(ancestor) get_n_synonym_mutations(ancestor))) 
  
  total_mutations <- node_mutations + ancestors_mutations
  return(total_mutations)
})

accumulated_mutations_descendants <- sapply(annotated_tree$node, function(n) {
  node_mutations <- get_n_synonym_mutations(n)
  descendants <- get_descendants(n)
  descendants_mutations <- sum(sapply(descendants, function(descendant) get_n_synonym_mutations(descendant))) 
  
  total_mutations <- node_mutations + descendants_mutations
  return(total_mutations)
})

annotated_tree$accumulated_mutations_root <- accumulated_mutations_root
annotated_tree$accumulated_mutations_descendants <- accumulated_mutations_descendants

