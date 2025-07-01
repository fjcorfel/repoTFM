library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(parallel)

CORES <- 10
YEARS_TO_ANALYZE <- 100

TOTAL_TREE_YEARS <- 4000

### FUNCTIONS -----------------------------------------------------------------

get_node_mutations <- function(node) {
  mutations <- unlist(result_tree[result_tree$node == node, "mutations"])
  if (length(mutations) == 0) return(character(0))
  
  return(mutations)
}

check_reversions <- function(node, mutation) {
  mutation_position <- as.numeric(str_extract(mutation, "\\d+"))
  descendants <- phangorn::Descendants(tree, node, "all")
  descendant_mutations <- unlist(sapply(descendants, get_node_mutations), use.names = FALSE)
  
  descendant_mutations <- na.omit(descendant_mutations)
  
  matches <- str_detect(descendant_mutations, paste0(mutation_position, "(?=\\D)"))
  return(any(matches))
}

check_mutations_in_wt <- function(node, mutation) {
  sister <- phangorn::Siblings(tree, node)
  sister_descendants <- unlist(phangorn::Descendants(tree, sister, "all"))
  sister_descendant_mutations <- unlist(sapply(sister_descendants, get_node_mutations))
  
  sister_descendant_mutations <- na.omit(sister_descendant_mutations)
  
  return(mutation %in% sister_descendant_mutations)
}

get_ancestors <- function(node) {
  ancestors <- phangorn::Ancestors(tree, node)
  return(ancestors)
}  

get_descendants <- function(node) {
  descendants <- phangorn::Descendants(tree, node, "all")
  return(descendants)
}  

get_n_synonym_mutations <- function(node) {
  mutation_idx <- which(result_tree$node == node)
  if (length(mutation_idx) > 0) {
    return(result_tree$n_synonym_mutations[mutation_idx])
  } else return(0)
}


### INPUTS --------------------------------------------------------------------

tree <- ape::as.phylo(treeio::read.beast("../../data/global/annotated_tree.nexus")) 

result_tree <- fread("../../data/global/annotated_tree_cleaned.csv") %>%
  as_tibble() %>%
  mutate(mutations = strsplit(mutations, "\\|"))


# print("Calculating accumulated mutations root")
# accumulated_mutations_root <- unlist(mclapply(result_tree$node, function(n) {
#   print(paste("Node root", n))
#   node_mutations <- get_n_synonym_mutations(n)
#   ancestors <- get_ancestors(n)
#   ancestors_mutations <- sum(unlist(sapply(ancestors, function(ancestor) get_n_synonym_mutations(ancestor)))) 
#   
#   total_mutations <- node_mutations + ancestors_mutations
#   return(total_mutations)
# }, mc.cores = 10, mc.preschedule = FALSE))
# save(accumulated_mutations_root, file = "../../data/global/global_accumulated_mutations_root.rda")
# print("ACUMULATED MUTATIONS 1 COMPLETED")
# 
# print("Calculating accumulated mutations descendants")
# accumulated_mutations_descendants <- unlist(mclapply(result_tree$node, function(n) {
#   print(paste("Node descendants", n))
#   node_mutations <- get_n_synonym_mutations(n)
#   descendants <- get_descendants(n)
#   descendants_mutations <- sum(unlist(sapply(descendants, function(descendant) get_n_synonym_mutations(descendant)))) 
#   
#   total_mutations <- node_mutations + descendants_mutations
#   return(total_mutations)
# }, mc.cores = 10, mc.preschedule = FALSE))
# save(accumulated_mutations_descendants, file = "../../data/global/global_accumulated_mutations_descendants.rda")
# print("ACUMULATED MUTATIONS 2 COMPLETED")
load("../../data/global/global_accumulated_mutations_descendants.rda")
load("../../data/global/global_accumulated_mutations_root.rda")

result_tree$accumulated_mutations_root <- accumulated_mutations_root
result_tree$accumulated_mutations_descendants <- accumulated_mutations_descendants


mutations_total_tree_years <- max(result_tree$accumulated_mutations_root)
mutations_years_to_analize <- round((YEARS_TO_ANALYZE * mutations_total_tree_years) / TOTAL_TREE_YEARS)

filtered_nodes <- result_tree %>%
  filter(accumulated_mutations_descendants <= mutations_years_to_analize)

nodes_branch_length <- fread("../../data/global/annotated_tree.csv") %>%
  select(node, branch.length) %>%
  as_tibble()

snp_count <- fread("../../data/global/SNP_count.csv") %>%
  filter(count >= 5)
mutations <- snp_count$mutation

snp_table <- fread("../../data/global/SNP_table_final.txt")


### MAIN PROCESSING -----------------------------------------------------------

process_mutation <- function(n_mutation) {
  
  print(paste("Processing mutation:", n_mutation))
  mutation <- mutations[n_mutation]
  
  # Filter nodes that have the mutation
  nodes_with_mutation <- filtered_nodes %>%
    rowwise() %>%
    filter(mutation %in% unlist(mutations)) %>%
    ungroup()
  
  # Verify that nodes are not empty
  if (nrow(nodes_with_mutation) == 0) {
    return(tibble(node=integer(), mutation=character(),
                  n_node_tips=integer(), n_sister_tips=integer()))
  }
  
  # Check if mutation and WT branches are pure
  pure_nodes <- nodes_with_mutation %>%
    rowwise() %>%
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
  if (nrow(pure_nodes) == 0) {
    return(tibble(node=integer(), mutation=character(),
                  n_node_tips=integer(), n_sister_tips=integer()))
  }
  
  final_nodes <- pure_nodes %>%
    mutate(
      mutation = mutation,
      RoHO = n_node_tips / n_sister_tips
    ) %>%
    select(node, mutation, n_node_tips, n_sister_tips, RoHO)
  
  return(final_nodes)
}


global_RoHO_nodes <- mclapply(seq_along(mutations), function(n_mutation) {
  process_mutation(n_mutation)
}, mc.cores = CORES, mc.preschedule = FALSE)

global_RoHO_nodes <- do.call(rbind, global_RoHO_nodes)

save(global_RoHO_nodes, file = paste0("../../data/global/global_RoHO_nodes_homoplasies_agefilter",YEARS_TO_ANALYZE,".rda"))

global_RoHO <- global_RoHO_nodes %>%
  na.omit() %>%
  group_by(mutation) %>%
  summarise(
    n_node_tips = sum(n_node_tips),
    n_sister_tips = sum(n_sister_tips),
    RoHO = n_node_tips / n_sister_tips
  ) %>%
  mutate(
    synonym = snp_table$Synonym[match(str_sub(mutation, 1, -2), snp_table$Position)],
    Rv_number = snp_table$Rv_number[match(str_sub(mutation, 1, -2), snp_table$Position)]
  ) %>%
  select(
    mutation, synonym, Rv_number, n_node_tips, n_sister_tips, RoHO
  ) %>%
  ungroup()

global_RoHO <- global_RoHO %>%
  mutate(synonym = ifelse(synonym %in% c("", "-"), NA, synonym))

save(global_RoHO, file = paste0("../../data/global/global_RoHO_homoplasies_agefilter",YEARS_TO_ANALYZE,".rda"))
fwrite(global_RoHO, file = paste0("../../data/global/global_RoHO_homoplasies_agefilter",YEARS_TO_ANALYZE,".csv"))
