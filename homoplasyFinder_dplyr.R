library(treeio)
library(phangorn)
library(dplyr)
library(data.table)
library(phytools)
library(parallel)


### LOAD INPUTS ###

print("Loading tree...")
tree <- ape::as.phylo(treeio::read.beast("../data/annotated_tree.nexus")) 

print("Loading SNP table mutations...")
snp_table <- data.table::fread("../data/SNP_table_noresis_mutations.txt")$Mutation

print("Loading ancestral mutations...")
load("../data/ancestral_result.rda")    # result_tree

# Convert to tibble for easier and more efficient dplyr manipulation
result_tree <- as_tibble(result_tree) %>%
  select(node, parent, label, ref_mutation_position)

print("Inputs loaded!")


### FUNCTIONS ###

# Get mutations from a given node
get_node_mutations <- function(node) {
  mutations <- unlist(result_tree[node, "ref_mutation_position"])
  if (is.null(mutations) || length(mutations) == 0) return(NULL)
  return(mutations)
}


# Check if mutation is in parent OR sister (mutation inherited)
# Return bool
mutation_in_sister_parent <- function(node, mutation) {
  
  if (is.null(node) || length(node) == 0) {
    return(FALSE)
  }
  
  sister_node <- phytools::getSisters(tree, node)
  sister_mutations <- get_node_mutations(sister_node)
  
  if (node == getRoot(tree)) {
    parent_mutations <- NULL
  } else {
    parent_node <- phytools::getParent(tree, node)
    parent_mutations <- get_node_mutations(parent_node)
  }
  
  return(mutation %in% c(sister_mutations, parent_mutations))
}

# Count MUT alleles
# Return int
count_mut_alleles <- function(tips, mutation) {
  tips <- unlist(tips)
  mut_alleles <- sum(sapply(tips, function(tip) mutation %in% get_node_mutations(tip)))
  return(mut_alleles)
}

count_wt_alleles <- function(tips, mut_alleles) {
  wt_alleles <- length(tips) - mut_alleles
  return(wt_alleles)
}

# Find if given SNP table position contains homoplasy
find_homoplasy <- function(n_position) {
  # Debugging
  if(n_position %% 1000 == 0) {
    print(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S")," - SNP Position: ", n_position))
  }
  
  # Get SNP mutation
  snp_position <- snp_table$Position[n_position]
  wt <- snp_table$WT[n_position]
  alt <- snp_table$ALT[n_position]
  snp_mutation <- paste0(wt, snp_position, alt)
  
  # Filter nodes that have the SNP mutation
  nodes_with_mutation <- result_tree %>%
    rowwise() %>%
    filter(snp_mutation %in% ref_mutation_position) %>%
    select(parent, node, label, ref_mutation_position) %>%
    ungroup()
  
  # Criteria for saving a node as homoplasy:
  # 1. Parent and sister don't have the mutation (mutation not inherited)
  # 2. At least 2 descendant tips of each allele (mutation and WT)
  
  # Filter nodes with the mutation in parent OR sister (mutation inherited)
  nodes_without_family_mutations <- nodes_with_mutation %>%
    rowwise() %>%
    mutate(
      mutation_in_family = mutation_in_sister_parent(node, snp_mutation)
    ) %>%
    filter(!mutation_in_family) %>%
    ungroup()
  
  # Filter nodes with 2 tips of each allele
  nodes_with_alleles <- nodes_without_family_mutations %>%
    rowwise() %>%
    mutate(
      node_tips = phangorn::Descendants(tree, node, "tips")
    ) %>%
    mutate(
      n_mut_alleles = count_mut_alleles(node_tips, snp_mutation)
    ) %>%
    mutate(
      n_wt_alleles = count_wt_alleles(node_tips, n_mut_alleles)
    ) %>% 
    filter(
      n_mut_alleles >= 2,
      n_wt_alleles >= 2
    ) %>%
    ungroup()
  
  homoplasy_nodes <- nodes_with_alleles %>%
    mutate(mutation = snp_mutation) %>%
    select(node, label, mutation, n_mut_alleles, n_wt_alleles)
    
  return(homoplasy_nodes)
}


### MAIN PROCESSING ###
start_time <- proc.time()
print("Starting parallel processing...")
n_cores <- detectCores() 

homoplasy_nodes <- mclapply(seq_along(snp_table), function(n_position) {
  find_homoplasy(n_position)
}, mc.cores = n_cores - 1)

# Combine df returned by each worker function into a single df
homoplasy_nodes <- do.call(rbind, homoplasy_nodes)

print("Processing complete.")

# Save the final result
print("Saving results...")
save(homoplasy_nodes, file = "../data/homoplasy_nodes.Rda")
#save(result_tree, file = "../data/ancestral_result_reversions.Rda")

end_time <- proc.time() - start_time

time_str <- paste("Runtime HomoplasyFinder:\n",
                  end_time["elapsed"], "secs")
writeLines(time_str, "runtime_HomoplasyFinder.txt")

print("HomoplasyFinder has finished!")



