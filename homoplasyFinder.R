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

print("Loading ancestral mutations...")
load("../data/ancestral_result.rda")    # result_tree

# Initialize column for tracking nodes containing reversions
result_tree$reversion <- FALSE


### HELPER FUNCTIONS ###

# Get mutations from a given node
get_node_mutations <- function(result_tree, node_number) {
  mutations <- unlist(result_tree[node_number, "ref_mutation_position"])
  if (is.null(mutations) || length(mutations) == 0) return(NULL)
  return(mutations)
}

# Check if mutation exists in BOTH sister and parent nodes
mutation_in_sister_parent <- function(tree, result_tree, node_number, mutation) {
  sister_node_number <- phytools::getSisters(tree, node_number)
  sister_mutations <- get_node_mutations(result_tree, sister_node_number)
  
  parent_node_number <- phytools::getParent(tree, node_number)
  parent_mutations <- get_node_mutations(result_tree, parent_node_number)
  
  return(mutation %in% c(sister_mutations, parent_mutations))
}

# Check if descendants acquired homoplasy or not (for checking posible reversions)
check_node_reversions <- function(tree, result_tree, node_number, mutation) {
  descendants <- unlist(phangorn::Descendants(tree, node_number, "all"))
  
  # Iterate through descendants to check for mutation loss
  for (descendant in descendants) {
    descendant_mutations <- get_node_mutations(result_tree, descendant)
    
    # If mutation is not present in this descendant, mark it as a reversion
    if (is.null(descendant_mutations) || !(mutation %in% descendant_mutations)){
      result_tree$reversion[descendant] <- TRUE
    }
  }
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
  print(paste0("Position: ", n_position))
  
  # Dataframe for saving homoplasy nodes that meet the criteria
  # 1. Parent and sister don't have the mutation (mutation not inherited)
  # 2. At least 2 descendant tips of each allele (mutation and WT)
  homoplasy_nodes <- data.frame(
    node_number = integer(),
    node_label = character(),
    node_mutation = character(),
    n_mut_alleles = integer(),
    n_wt_alleles = integer(),
    stringsAsFactors = FALSE
  )
  
  # Get SNP mutation
  snp_position <- snp_table$Position[n_position]
  wt <- snp_table$WT[n_position]
  alt <- snp_table$ALT[n_position]
  snp_mutation <- paste0(wt, snp_position, alt)
  

  for (n_node in seq_along(result_tree$node)) {
    # Check if node has the snp_mutation
    node_mutations <- get_node_mutations(result_tree, n_node)
    if (is.null(node_mutations) || !snp_mutation %in% node_mutations) next
    
    # Check if sister and parent have the snp_mutation
    if (mutation_in_sister_parent(tree, result_tree, n_node, snp_mutation)) next
    
    # Check alleles in tips
    tip_allele_counts <- count_tip_alleles(tree, result_tree, n_node, snp_mutation)
    if (is.null(tip_allele_counts)) next
    
    # Check if descendant nodes mantain mutation or not (reversion)
    check_node_reversions(tree, result_tree, n_node, snp_mutation)
    
    
    # If node meets all the criteria -> save node as row in df
    homoplasy_nodes <- rbind(
      homoplasy_nodes,
      data.frame(
        node_number = n_node,
        node_label = result_tree[n_node, "label"],
        node_mutation = snp_mutation,
        n_mut_alleles = tip_allele_counts$mut_alleles,
        n_wt_alleles = tip_allele_counts$wt_alleles,
        stringsAsFactors = FALSE
      )
    )
  }
  
  return(homoplasy_nodes)
}



### MAIN PROCESSING ###

print("Starting parallel processing...")
n_cores <- detectCores() 

homoplasy_nodes <- mclapply(seq_along(snp_table$Position), function(n_position) {
  find_homoplasy(n_position, snp_table, tree, result_tree)
}, mc.cores = 8)

# Combine df returned by each worker function into a single df
homoplasy_nodes <- do.call(rbind, homoplasy_nodes)

# Save the final result
save(homoplasy_nodes, file = "homoplasy_nodes.Rda")

print("Processing complete.")
