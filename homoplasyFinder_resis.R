library(treeio)
library(phangorn)
library(dplyr)
library(data.table)
library(phytools)
library(parallel)
library(stringr)


### FILTER SNPS ###

# Load SNP count
load("../data/SNP_count_resis.rda")
# Filter SNPs that appear more than once
snp_count_filtered <- snp_count[snp_count >= 5]

# Load SNP table and select mutations
snp_table_mutations <- data.table::fread("../data/SNP_table_resis.txt") %>%
  select(Position, WT, ALT) %>%
  mutate(Mutation = paste0(WT, Position, ALT)) %>%
  select(Mutation)

filtered_mutations <- snp_table_mutations$Mutation[snp_table_mutations$Mutation %in% names(snp_count_filtered)]
save(filtered_mutations, file = "../data/filtered_mutations_resis.rda")


### LOAD INPUTS ###

# Load tree -> phylo object
print("Loading tree...")
tree <- ape::as.phylo(treeio::read.beast("../data/annotated_tree_resis.nexus")) 

# Load SNP table mutations -> vector
# Mutations previously filtered for > 1 occurrence (possible homoplasy)
print("Loading SNP mutations...")
load("../data/filtered_mutations_resis.rda")

# Load tree nodes and their associated mutations -> tibble
print("Loading ancestral nodes and mutations...")
load("../data/ancestral_result_resis.rda")    # result_tree
# Convert to tibble for more efficient dplyr manipulation
result_tree_resis <- as_tibble(result_tree_resis) %>%
  select(node, parent, label, ref_mutation_position)

print("Inputs loaded.")


### FUNCTIONS ###

# Get mutations from a given node
get_node_mutations <- function(node) {
  mutations <- unlist(result_tree_resis[node, "ref_mutation_position"])
  if (is.null(mutations) || length(mutations) == 0) return(NULL)
  return(mutations)
}

# Get mut alleles from a given node (check node tips )
get_mut_alleles <- function(node) {
  tips <- unlist(phangorn::Descendants(tree, node, "tips"))
  
  return(length(tips))
}

# Get wt alleles from a given node (check sister tips)
get_wt_alleles <- function(node) {
  sister <- phytools::getSisters(tree, node)
  sister_tips <- unlist(phangorn::Descendants(tree, sister, "tips"))
  
  return(length(sister_tips))
}

check_reversions <- function(node, mutation) {
  mutation_position <- as.numeric(str_extract(mutation, "\\d+"))
  descendants <- phangorn::Descendants(tree, node, "all")
  descendant_mutations <- unlist(sapply(descendants, get_node_mutations), use.names = FALSE)
  
  matches <- str_detect(descendant_mutations, paste0("(?<=\\D)", mutation_position, "(?=\\D)"))
  return(any(matches))
}

check_mutations_in_wt <- function(node, mutation) {
  sister <- phytools::getSisters(tree, node)
  sister_descendants <- unlist(phangorn::Descendants(tree, sister, "all"))
  sister_descendant_mutations <- unlist(sapply(sister_descendants, get_node_mutations))
  
  return(mutation %in% sister_descendant_mutations)
}


# Find if given SNP table position contains homoplasy
find_homoplasy <- function(n_position) {
  
  # Get SNP mutation
  snp_mutation <- filtered_mutations[n_position]
  
  # Debugging
  if(n_position %% 1 == 0) {
    print(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                 " - Position: ", n_position,
                 " - Mutation: ", snp_mutation))
  }
  
  # Filter nodes that have the SNP mutation
  nodes_with_mutation <- result_tree_resis %>%
    rowwise() %>%
    filter(snp_mutation %in% unlist(ref_mutation_position)) %>%
    select(parent, node, label, ref_mutation_position) %>%
    ungroup()
  
  
  # Criteria for saving a node as homoplasy:
  # WT and mutation branches are pure (no appearing of the mutation or reversions)
  # At least 2 descendant tips of each allele (mutation and WT)
  
  # Filter nodes with 2 tips of each allele
  nodes_with_alleles <- nodes_with_mutation %>%
    rowwise() %>%
    mutate(
      n_mut_alleles = get_mut_alleles(node)
    ) %>%
    mutate(
      n_wt_alleles = get_wt_alleles(node)
    ) %>%
    filter(
      n_mut_alleles >= 2,
      n_wt_alleles >= 2
    ) %>%
    ungroup()
  
  # Verify that nodes are not empty before continue
  if (nrow(nodes_with_alleles) == 0) {
    return(tibble(node=integer(), label = character(), mutation = character(),
                  n_mut_alleles = integer(), n_wt_alleles = integer()))
  }
  
  # Check if mutation and WT branches are pure
  filtered_alleles <- nodes_with_alleles %>%
    rowwise() %>%
    mutate(
      has_reversion = check_reversions(node, snp_mutation),
      has_mutation_in_wt = check_mutations_in_wt(node, snp_mutation)
    ) %>%
    filter(
      !has_reversion, 
      !has_mutation_in_wt
    ) %>%
    ungroup()
  
  # Verify that nodes are not empty before continue
  if (nrow(filtered_alleles) == 0) {
    return(tibble(node=integer(), label = character(), mutation = character(),
                  n_mut_alleles = integer(), n_wt_alleles = integer()))
  }
  
  homoplasy_nodes <- filtered_alleles %>%
    mutate(
      mutation = snp_mutation,
      RoHO = n_mut_alleles / n_wt_alleles
    ) %>%
    select(node, label, mutation, n_mut_alleles, n_wt_alleles, RoHO)
  
  return(homoplasy_nodes)
}


### MAIN PROCESSING ###

print("Starting parallel processing...")
n_cores <- detectCores() 

homoplasy_nodes <- mclapply(seq_along(filtered_mutations), function(n_position) {
  find_homoplasy(n_position)
}, mc.cores = 14)

# Combine df returned by each worker function into a single df
homoplasy_nodes_resis <- do.call(rbind, homoplasy_nodes)

print("Processing complete.")
# Remove result tree from environment for saving memory


### SAVING RESULTS ### 

# Save the final result
print("Saving results...")
save(homoplasy_nodes_resis, file = "../data/homoplasy_nodes_resis.rda")

print("HomoplasyFinder has finished!")
