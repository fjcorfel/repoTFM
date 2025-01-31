library(treeio)
library(phangorn)
library(dplyr)
library(data.table)
library(phytools)
library(parallel)
library(stringr)


### LOAD INPUTS ###

# Load tree -> phylo object
print("Loading tree...")
tree <- ape::as.phylo(treeio::read.beast("../data/annotated_tree.nexus")) 

# Load SNP table mutations -> vector
# Mutations previously filtered for > 1 occurrence (possible homoplasy)
print("Loading SNP mutations...")
load("../data/filtered_mutations.rda")

# Load tree nodes and their associated mutations -> tibble
print("Loading ancestral nodes and mutations...")
load("../data/ancestral_result.rda")    # result_tree
# Convert to tibble for more efficient dplyr manipulation
result_tree <- as_tibble(result_tree) %>%
  select(node, parent, label, ref_mutation_position)

print("Inputs loaded.")


### FUNCTIONS ###

# Get mutations from a given node
get_node_mutations <- function(node) {
  mutations <- unlist(result_tree[node, "ref_mutation_position"])
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
  if(n_position %% 1000 == 0) {
    print(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                 " - Position: ", n_position,
                 " - Mutation: ", snp_mutation))
  }
  
  # Filter nodes that have the SNP mutation
  nodes_with_mutation <- result_tree %>%
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
homoplasy_nodes <- do.call(rbind, homoplasy_nodes)

print("Processing complete.")
# Remove result tree from environment for saving memory
rm(result_tree)

# Save the final result
print("Saving results...")
save(homoplasy_nodes, file = "../data/homoplasy_nodes.rda")

print("Annotating mutations...")
# Read SNP table (complete version)
snp_table <- as_tibble(fread("../data/SNP_table_noresis.txt"))  
snp_table <- snp_table %>%
  select(Position, WT, ALT, Synonym, Rv_number) %>%
  mutate(Mutation = paste0(WT, Position, ALT)) %>%
  select(Mutation, Synonym, Rv_number)

# Add gene names and Rv number 
homoplasy_nodes_annotated <- homoplasy_nodes %>%
  mutate(
    synonym = snp_table$Synonym[match(mutation, snp_table$Mutation)],
    Rv_number = snp_table$Rv_number[match(mutation, snp_table$Mutation)]
  ) %>%
  select(
    node,
    label,
    mutation,
    synonym,
    Rv_number,
    n_mut_alleles,
    n_wt_alleles,
    RoHO)

# Save results
save(homoplasy_nodes_annotated, file = "../data/homoplasy_nodes_annotated.rda")

# Group by mutation
homoplasy_nodes_annotated_byMutation <- homoplasy_nodes_annotated %>%
  group_by(mutation) %>%
  summarise(
    n_mut_alleles = sum(n_mut_alleles),
    n_wt_alleles = sum(n_wt_alleles),
    RoHO = n_mut_alleles / n_wt_alleles
  ) %>%
  mutate(
    synonym = snp_table$Synonym[match(mutation, snp_table$Mutation)],
    Rv_number = snp_table$Rv_number[match(mutation, snp_table$Mutation)]
  ) %>%
  select(
    mutation,
    synonym,
    Rv_number,
    n_mut_alleles,
    n_wt_alleles,
    RoHO
  ) %>%
  ungroup()

save(homoplasy_nodes_annotated_byMutation, file = "../data/homoplasy_mutations.rda")

# Group by gene
homoplasy_nodes_annotated_byGene <- homoplasy_nodes_annotated %>%
  group_by(Rv_number) %>%
  summarise(
    n_mut_alleles = sum(n_mut_alleles),
    n_wt_alleles = sum(n_wt_alleles),
    RoHO = n_mut_alleles / n_wt_alleles
  ) %>%
  mutate(
    synonym = snp_table$Synonym[match(Rv_number, snp_table$Rv_number)]
  ) %>%
  select(
    Rv_number,
    synonym,
    n_mut_alleles,
    n_wt_alleles,
    RoHO
  ) %>%
  ungroup()

save(homoplasy_nodes_annotated_byGene, file = "../data/homoplasy_genes.rda")

print("HomoplasyFinder has finished!")
