library(treeio)
library(phangorn)
library(dplyr)
library(data.table)
library(phytools)

### Load inputs ###

# Tree
print("Loading tree...")
tree <- ape::as.phylo(treeio::read.beast("../data/annotated_tree.nexus")) 
# SNP table
print("Loading SNP table...")
snp_table <- data.table::fread("../data/SNP_table_noresis.txt") %>%
  select(Position, WT, ALT)
print("Loading mutations table (ancestral recosntruction)...")
# result_tree (mutations per node table)
load("../data/ancestral_result.rda")


### Processing ###

# List for saving nodes associated to each SNP mutation
node_mutations_list <- list()

# Iterate for each position of SNP table
for (n_position in seq_along(snp_table$Position)) {
  
  # Extract mutation
  position <- snp_table$Position[n_position]
  wt <- snp_table$WT[n_position]
  alt <- snp_table$ALT[n_position]
  snp_mutation <- paste0(wt, position, alt)
  
  # Iterate for each node of ancestral result
  for (n_node in seq_along(result_tree$node)) {
    
    # Get mutations of the node
    ref_mutations <- unlist(result_tree[n_node, "ref_mutation_position"])
    if (is.null(ref_mutations) || length(ref_mutations) == 0) next
    
    # Iterate over mutations of the node
    for (node_mutation in ref_mutations) {
      
      # If node has the SNP mutation, we save it
      if (node_mutation == snp_mutation) {
        node_label <- result_tree[n_node, ]$label
        node_info <- list(
          node_number = n_node,
          node_label = node_label,
          node_mutation = snp_mutation
        )
        node_mutations_list <- append(node_mutations_list, list(node_info))
        
      }
    }
  }
}

### Filtering ###

# List of nodes that meet the filtering criteria
homoplasy_nodes <- list()

# Iterate over saved nodes for filtering
for (n_node in seq_along(node_mutations_list)) {
  
  node_info <- node_mutations_list[[n_node]]
  
  # Check if sister or parent has mutation
  sister_node_number <- phytools::getSisters(tree, node_info$node_number)
  sister_node_mutations <- unlist(result_tree[sister_node_number, "ref_mutation_position"])
  
  parent_node_number <- phytools::getParent(tree, node_info$node_number)
  parent_node_mutations <- unlist(result_tree[parent_node_number, "ref_mutation_position"])
  
  if (mutation %in% sister_node_mutations || mutation %in% parent_node_mutations) next # Skip this node
  
  # Check if node has at least 2 tips for each allele
  tips <- unlist(phangorn::Descendants(tree, node_info$node_number, "tips"))
  
  if (length(tips) < 4) next # Skip this node
  
  mut_alleles <- 0
  wt_alleles <- 0
  for (tip in tips) {
    tip_mutations <- unlist(result_tree[tip, "ref_mutation_position"])
    
    if (mutation %in% tip_mutations) {
      mut_alleles <- mut_alleles + 1
    } else {
      wt_alleles <- wt_alleles + 1
    }
  }
  
  # We save the node and its tips info if it meets the criteria
  node_info[["mut_alleles"]] <- mut_alleles
  node_info[["wt_alleles"]] <- wt_alleles
  
  # NODOS CON MUTACIONES QUE SUPEREN EL FILTRADO SE MARCARÁN COMO HOMOPLASIAS
  homoplasy_nodes <- append(homoplasy_nodes, list(node_info))
}


  
