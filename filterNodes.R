library(data.table)
library(dplyr)
library(tidyr)
library(phangorn)
library(parallel)
library(stringr)


DATASET <- "vietnam"
DATASET_SIZE <- switch (DATASET,
  "malawi" = 1840,
  "vietnam" = 1504
)


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

calculate_n_sibling_tips <- function(node) {
  siblings <- unlist(phangorn::Siblings(tree, node))
  sibling_tips <- unlist(phangorn::Descendants(tree, siblings, "tips"))
  return(length(sibling_tips))
}

get_node_branch_length <- function(node) {
  node_branch_length <- nodes_branch_length[[node, "branch.length"]]
}

test_branch_lengths <- function(node) {
  node_tips <- unlist(phangorn::Descendants(tree, node, "tips"))
  siblings <- unlist(phangorn::Siblings(tree, node))
  sibling_tips <- unlist(phangorn::Descendants(tree, siblings, "tips"))
  
  mutant_branch_lengths <- sapply(node_tips, get_node_branch_length)
  wt_branch_lengths <- sapply(sibling_tips, get_node_branch_length)
  
  # This block handles constant data error in t.test
  # This error happens when there is not branch length variation
  results <- tryCatch({
    t.test(mutant_branch_lengths, wt_branch_lengths, alternative = "less")
  }, error = function(e) {
    return(NA)
  })
  
  if (inherits(results, "htest")) {
    return(results$p.value)
  } else {
    return(NA)
  }
}

### INPUTS --------------------------------------------------------------------

# Load tree (nwk)
tree <- ape::as.phylo(treeio::read.beast(paste0("../data/", DATASET, "/", "annotated_tree.nexus"))) 


# Load annotated_tree_cleaned (mutations per node)

annotated_tree <- fread(paste0("../data/", DATASET, "/", "annotated_tree_cleaned.csv")) %>%
  as_tibble() %>%
  mutate(mutations = strsplit(mutations, '\\|')) %>%
  rowwise() %>%
  mutate(n_tips = calculate_n_tips(node),
         n_sibling_tips = calculate_n_sibling_tips(node)) %>%
  filter((n_tips >= 2 & n_tips <= (DATASET_SIZE / 10)) & 
           (n_sibling_tips >= 2 & n_sibling_tips <= (DATASET_SIZE / 10))) %>%
  ungroup()

# Load annotated_tree (branch length per node)

nodes_branch_length <- fread(paste0("../data/", DATASET, "/", "annotated_tree.csv")) %>%
  select(node, branch.length) %>%
  as_tibble()

# Select mutations with more than 1 appearence

snp_count <- fread(paste0("../data/", DATASET, "/", "SNP_count.csv")) %>%
  filter(count > 1)
mutations <- snp_count$mutation


### MAIN PROCESSING -----------------------------------------------------------

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
    ungroup()
  
  # Verify that nodes are not empty
  if (nrow(nodes) == 0) {
    return(tibble(node=integer(), mutation=character(),branch_length=character()))
  }
  
  nodes <- nodes %>%
    rowwise() %>%
    mutate(
      has_reversion = check_reversions(node, mutation),
      has_mutation_in_wt = check_mutations_in_wt(node, mutation)
    ) %>%
    ungroup() %>%
    filter(
      !has_reversion,
      !has_mutation_in_wt
    ) 
  
  # Verify that nodes are not empty
  if (nrow(nodes) == 0) {
    return(tibble(node=integer(), mutation=character(),branch_length=character()))
  }
  
  # Hacer aquí el cálculo
  # Sacar las longitudes de tips mutantes y tips wt
  # Hacer t.test entre ambos vectores -> guardar pvalor
  filtered_nodes <- nodes %>%
    mutate(mutation = mutation) %>%
    rowwise() %>%
    mutate(ttest_pvalue = test_branch_lengths(node)) %>%
    ungroup() %>%
    mutate(ttest_adj_pvalue_Bonf = p.adjust(ttest_pvalue, method = "bonferroni"),
           ttest_adj_pvalue_BH = p.adjust(ttest_pvalue, method = "BH")) %>%
    select(node, mutation, ttest_pvalue, ttest_adj_pvalue_Bonf, ttest_adj_pvalue_BH)
  
  return(filtered_nodes)
  
}

final_nodes <- mclapply(seq_along(mutations), function(n_mutation) {
  tryCatch({
    process_mutation(n_mutation)
  }, error = function(e) {
    write(paste("Error en mutación:", mutations[n_mutation], "->", conditionMessage(e)),
          file = "error_log.txt", append = TRUE)
    return(NULL)
  })
  
}, mc.cores = 10, mc.preschedule = FALSE)

final_nodes <- do.call(rbind, final_nodes)

fwrite(final_nodes, file = paste0("../data/", DATASET, "/", "final_nodes_", DATASET, ".csv"))

# Group by mutation
final_mutations <- final_nodes %>%
  drop_na() %>%
  group_by(mutation) %>%
  summarise(adj_pvalues = list(ttest_adj_pvalue_BH),
            n = n(),
            significant_ratio = sum(ttest_adj_pvalue_BH <= 0.05) / n()) %>%
  ungroup()
  
fwrite(final_mutations, file = paste0("../data/", DATASET, "/", "final_mutations_", DATASET, ".csv"))
