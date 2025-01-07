library(treeio)
library(phangorn)
library(dplyr)
library(data.table)
library(phytools)

### Load inputs ###
print("--- LOADING INPUTS ---")
# Tree
print("Loading tree...")
tree <- ape::as.phylo(treeio::read.beast("../data/annotated_tree.nexus")) 
# SNP table
print("Loading SNP table")
snp_table <- data.table::fread("../data/SNP_table_noresis.txt") %>%
  select(Position, WT, ALT)
print("Loading mutations table (ancestral recosntruction)...")
# result_tree (mutations per node table)
load("../data/ancestral_result.rda")


### Processing nodes ###
print("--- PROCESSING NODES ---")

# List for saving results
nodes_with_mutations_info <- list()

# Iterate for each position of SNP table
for (n_position in seq_along(snp_table$Position)) {
  print("--------------------------------")
  print(paste("SNP POSITION:", n_position))
  
  # Extract mutation
  position <- snp_table$Position[n_position]
  wt <- snp_table$WT[n_position]
  alt <- snp_table$ALT[n_position]
  snp_mutation <- paste0(wt, position, alt)
  print(paste("MUTATION:", snp_mutation))
  
  
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
        nodes_with_mutations_info <- append(nodes_with_mutations_info, list(node_info))
        
      }
    }
  }
}

# List of nodes that pass the filtering criteria
homoplasy_nodes <- list()

# Iterate over saved nodes for filtering
for (n_node in seq_along(nodes_with_mutations_info)) {
  
  node_info <- nodes_with_mutations_info[[n_node]]
  
  # Check if descendants have the mutation position (check reversions)
  node_descendants <- phytools::getDescendants(tree, node_info$node_number)
  mutation <- node_info$node_mutation
  
  reversion_found <- FALSE
  for (descendant in node_descendants) {
    descendant_mutations <- unlist(result_tree[descendant, "ref_mutation_position"])
    if (mutation %in% descendant_mutations) {
      reversion_found <- TRUE
      # AÑADIR REVERSIÓN A LA TABLA (añadir columna en SNP table)
      break
    }
  }
  
  if (reversion_found) next # Skip this node
  
  
  # Check if sister or parent has mutation
  sister_node_number <- phytools::getSisters(tree, node_info$node_number)
  sister_node_mutations <- unlist(result_tree[sister_node_number, "ref_mutation_position"])
  
  parent_node_number <- phytools::getParent(tree, node_info$node_number)
  parent_node_mutations <- unlist(result_tree[parent_node_number, "ref_mutation_position"])
  
  if (mutation %in% sister_node_mutations || mutation %in% parent_node_mutations) next # Skip this node
  
  # Check if node has at least 2 tips for each allele
  tips <- unlist(phangorn::Descendants(tree, node_info$node_number, "tips"))
  
  if (length(tips) < 4) next # Skip this node
  
  
  #? Como se puede comprobar el alelo de cada tip ¿?
  
  # NODOS CON MUTACIONES QUE SUPEREN EL FILTRADO SE MARCARÁN COMO HOMOPLASIAS
  homoplasy_nodes <- append(homoplasy_nodes, list(node_info))
}


  
