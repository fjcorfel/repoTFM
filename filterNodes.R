library(data.table)
library(dplyr)
library(phangorn)


### INPUTS --------------------------------------------------------------------

DATASET_SIZE <- 1840

# Load tree (nwk)
tree <- ape::as.phylo(treeio::read.beast("../data/annotated_tree.nexus")) 

# Load annotated_tree_cleaned (mutations per node)
annotated_tree <- fread('../data/annotated_tree_cleaned.csv') %>%
  as_tibble() %>%
  mutate(mutations = strsplit(mutations, '\\|')) %>%
  rowwise() %>%
  mutate(n_tips = calculate_n_tips(node)) %>%
  filter(n_tips >= 2 & n_tips <= (DATASET_SIZE / 10)) %>%
  ungroup()
  

# Load annotated_tree (branch length per node)
nodes_branch_length <- fread('../data/annotated_tree.csv') %>%
  select(node, branch.length) %>%
  as_tibble()

# Select mutations with more than 1 appearence
snp_count <- fread('../data/SNP_count.csv') %>%
  filter(count > 1)
mutations <- snp_count$mutation


### FUNCTIONS -----------------------------------------------------------------

get_node_mutations <- function(node) {
  if (!node %in% annotated_tree$node) return(character(0))
  
  node_mutations <- unlist(annotated_tree[[which(annotated_tree$node == node), "mutations"]])
  return(node_mutations)
}


check_reversions <- function(node, mutation) {
  mutation_position <- as.numeric(str_extract(mutation, "\\d+"))
  
  descendants <- phangorn::Descendants(tree, node, "all")
  descendant_mutations <- unlist(lapply(descendants, get_node_mutations), use.names = FALSE)
  
  matches <- str_detect(descendant_mutations, paste0(mutation_position, "(?=\\D)"))
  return(any(matches))
}

check_mutations_in_wt <- function(node, mutation) {
  siblings <- phangorn::Siblings(tree, node)
  sibling_descendants <- unlist(phangorn::Descendants(tree, siblings, "all"))
  sibling_descendant_mutations <- unlist(lapply(sibling_descendants, get_node_mutations), use.names = FALSE)
  
  return(mutation %in% sibling_descendant_mutations)
}

calculate_n_tips <- function(node) {
  tips <- unlist(phangorn::Descendants(tree, node, "tips"))
  return(length(tips))
}



# Iterate over mutations
process_mutation <- function(n_mutation) {
  
  print(paste('Processing mutation:', n_mutation))
  
  mutation <- mutations[n_mutation]
  
  # Filter nodes that have correct tip number size
  # Filter nodes that have the mutation
  # Check if mutation and WT branches are pure
  nodes <- annotated_tree %>%
    rowwise() %>%
    filter(mutation %in% unlist(mutations)) %>%
    mutate(
      has_reversion = check_reversions(node, mutation),
      has_mutation_in_wt = check_mutations_in_wt(node, mutation)
    ) %>%
    filter(
      !has_reversion,
      !has_mutation_in_wt
    ) %>%
    ungroup()
  
  # Verify that nodes are not empty
  if (nrow(nodes) == 0) {
    return(tibble(node=integer(), mutation=character(),branch_length=character()))
  }
  
  final_nodes <- nodes %>%
    mutate(mutation = mutation) %>%
    select(node, mutation, branch.length)
  
  #? Lo que nosotros queremos es la longitud de rama de las tips
  
  return(final_nodes)
  
}

filtered_nodes <- mclapply(seq_along(mutations), function(n_mutation) {
  process_mutation(n_mutation)
}, mc.cores = 8, mc.preschedule = FALSE)

filtered_nodes <- do.call(rbind, filtered_nodes)

#?
# Get node tips and their branch length
# Get sister tips and their branch length
