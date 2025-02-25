library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(ggplot2)
library(parallel)


result_tree <- fread("../data/ancestral_result.csv")
result_tree <- result_tree %>%
  mutate(mutations = ifelse(mutations == "", NA, strsplit(mutations, "\\|"))) %>%
  as_tibble()

snp_table <- as_tibble(fread("../data/SNP_table_final_redundant.txt"))
tree <- ape::as.phylo(treeio::read.beast("../data/annotated_tree_resis.nexus")) 
load("../data/n_node_tips.rda")
load("../data/n_sister_tips.rda")


get_node_mutations <- function(node) {
  mutations <- unlist(result_tree[[node, "mutations"]])
  if (length(mutations) == 0) return(NA)

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


# Only with nodes with 4 < n_tips (samples) < 100
filtered_nodes <- result_tree %>% 
  mutate(n_node_tips = n_node_tips,
         n_sister_tips = n_sister_tips) %>%
  filter(n_node_tips > 4 & n_node_tips < 100)
 

# Extract all the mutations to process
mutations <- unique(unlist(filtered_nodes$mutations))

process_mutation <- function(n_mutation) {
  
  print(paste("Processing mutation:", n_mutation))
  mutation <- mutations[n_mutation]
  
  # Filter nodes that have the mutation
  nodes_with_mutation <- filtered_nodes %>%
    rowwise() %>%
    filter(mutation %in% unlist(mutations)) %>%
    ungroup()
  
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
}, mc.cores = 8, mc.preschedule = FALSE)

global_RoHO_nodes <- do.call(rbind, global_RoHO_nodes)

save(global_RoHO_nodes, file = "../data/global_RoHO_nodes.rda")

# global_RoHO <- tree_tibble %>%
#   mutate(n_node_tips = n_node_tips,
#          n_sister_tips = n_sister_tips) %>%
#   filter(n_node_tips > 4 & n_node_tips < 100) %>%
#   mutate(RoHO = n_node_tips / n_sister_tips) %>%
#   mutate(mutations = ifelse(
#     result_tree_noresis$ref_mutation_position[match(node, result_tree_noresis$node)] != "",
#     result_tree_noresis$ref_mutation_position[match(node, result_tree_noresis$node)],
#     result_tree_resis$ref_mutation_position[match(node, result_tree_resis$node)]
#   )) %>%
#   filter(mutations != "") %>%
#   separate_rows(mutations, sep = "\\|") %>%
#   rename(mutation = mutations) %>%
#   group_by(mutation) %>%
#   summarise(
#     n_node_tips = sum(n_node_tips),
#     n_sister_tips = sum(n_sister_tips),
#     RoHO = n_node_tips / n_sister_tips
#   ) %>%
#   mutate(
#     synonym = snp_table$Synonym[match(str_sub(mutation, 1, -2), snp_table$Position)],
#     Rv_number = snp_table$Rv_number[match(str_sub(mutation, 1, -2), snp_table$Position)]
#   ) %>%
#   select(
#     mutation, synonym, Rv_number, n_node_tips, n_sister_tips, RoHO
#   )
# 
# save(global_RoHO, file = "../data/global_RoHO.rda")
# fwrite(global_RoHO, file = "../data/global_RoHO.csv")
# 
# histogram <- global_RoHO %>%
#   ggplot(aes(x = log2(RoHO))) +
#   geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.4) +
#   ggtitle("Distribution of RoHO values")
# 
# histogram
