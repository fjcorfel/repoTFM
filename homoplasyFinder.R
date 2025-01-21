library(treeio)
library(phangorn)
library(dplyr)
library(data.table)
library(phytools)
library(parallel)


### LOAD INPUTS ###

print("Loading tree...")
tree <- ape::as.phylo(treeio::read.beast("../data/annotated_tree.nexus")) 

print("Loading SNP table...")
snp_table <- data.table::fread("../data/SNP_table_noresis.txt") %>%
  select(Position, WT, ALT)

print("Loading mutations table (ancestral recosntruction)...")
load("../data/ancestral_result.rda")



### HELPER FUNCTIONS ###

# Get mutations from a given node
get_node_mutations <- function(result_tree, node_number) {
  mutations <- unlist(result_tree[node_number, "ref_mutation_position"])
  if (is.null(mutations) || length(mutations) == 0) return(NULL)
  return(mutations)
}

# Check if mutation exists in sister AND parent nodes
mutation_in_sister_parent <- function(tree, result_tree, node_number, mutation) {
  sister_node_number <- phytools::getSisters(tree, node_number)
  sister_mutations <- get_node_mutations(result_tree, sister_node_number)
  
  parent_node_number <- phytools::getParent(tree, node_number)
  parent_mutations <- get_node_mutations(result_tree, parent_node_number)
  
  return(mutation %in% c(sister_mutations, parent_mutations))
}

#? Check if descendants acquired homoplasy (check reversions)
mutation_in_descendants <- function(tree, result_tree, node_number, mutation) {
  descendants <- unlist(phangorn::Descendants(tree, node_number, "all"))
  if (is.null(descendants) || length(descendants) == 0) return (FALSE)
  
  any(sapply(descendants, function(descendant) {
    mutations <- unlist(result_tree[descendant, "ref_mutation_position"])
    !is.null(mutations) && mutation %in% mutations
  }))
  
  #! TODO: habría que añadir FLAG a result_tree
}

# Count alleles for mutation and WT in tips
count_tip_alleles <- function(tree, result_tree, node_number, mutation) {
  tips <- unlist(phangorn::Descendants(tree, node_number, "tips"))
  
  # We need at least 2 tips with each allele
  if (length(tips) < 4) return(NULL)
  
  mut_alleles <- sum(sapply(tips, function(tip) mutation %in% get_node_mutations(result_tree, tip)))
  wt_alleles <- length(tips) - mut_alleles
  
  return(list(mut_alleles = mut_alleles, wt_alleles = wt_alleles))
}

# Process a SNP position from SNP table across all nodes
# Worker function for parallelization
find_homoplasy <- function(n_position, snp_table, tree, result_tree) {
  
  # List for saving homoplasy nodes that meet the criteria
  # 1. No descendant nodes that acquired homoplasy
  # 2. At least 2 descendant tips of each allele
  homoplasy_nodes <- list()
  
  # Get SNP mutation
  snp_position <- snp_table$Position[n_position]
  wt <- snp_table$WT[n_position]
  alt <- snp_table$ALT[n_position]
  snp_mutation <- paste0(wt, snp_position, alt)
  

  for (n_node in seq_along(result_tree$node)) {
    # Check if node has the snp_mutation
    node_mutations <- get_node_mutations(result_tree, n_node)
    if (is.null(node_mutations) || !snp_mutation %in% node_mutations) next
    
    # Check sister and parent nodes
    if (mutation_in_sister_parent(tree, result_tree, n_node, snp_mutation)) next
    
    # Check alleles in tips
    tip_allele_counts <- count_tip_alleles(tree, result_tree, n_node, snp_mutation)
    if (is.null(tip_allele_counts)) next
    
    # If node meets all the criteria -> save node
    node_info <- list(
      node_number = n_node,
      node_label = result_tree[n_node, "label"],
      node_mutation = snp_mutation,
      mut_alleles = tip_allele_counts$mut_alleles,
      wt_alleles = tip_allele_counts$wt_alleles
    )
    homoplasy_nodes <- append(homoplasy_nodes, list(node_info))
  }
  
  return(homoplasy_nodes)
}



### MAIN PROCESSING ###

print("Starting parallel processing...")
n_cores <- detectCores() 

homoplasy_nodes <- mclapply(seq_along(snp_table$Position), function(n_position) {
  
  find_homoplasy(n_position, snp_table, tree, result_tree)
  
}, mc.cores = n_cores)

# Combine lists returned by each worker function into a single list
homoplasy_nodes <- do.call(c, homoplasy_nodes)


print("Processing complete.")