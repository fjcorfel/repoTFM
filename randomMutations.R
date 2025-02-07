library(dplyr)
library(parallel)
library(pbapply)   # Progress bar for parallel loops
library(phangorn)
library(phytools)
library(rvmethod)
library(ggplot2)

### PRE PROCESSING ############################################################

# # Set up parallel processing
# num_cores <- detectCores() - 1  # Use max cores minus 1
# cl <- makeCluster(num_cores)
# clusterExport(cl, c("tree", "tree_tibble"))  # Export necessary objects to workers
# 
# # Compute Descendants in Parallel
# n_node_tips <- pbsapply(tree_tibble$node, function(n) length(unlist(phangorn::Descendants(tree, n, "tips"))), cl = cl)
# 
# # Compute Sisters in Parallel
# sisters <- pbsapply(tree_tibble$node, function(n) phytools::getSisters(tree, n), cl = cl)
# 
# # Compute Sister Descendants in Parallel
# n_sister_tips <- pbsapply(sisters, function(s) {
#   if (is.null(s)) return(0)  # Handle nodes with no sister
#   length(unlist(phangorn::Descendants(tree, s, "tips")))
# }, cl = cl)
# 
# # Stop parallel cluster
# stopCluster(cl)
# 

### FUNCTIONS #################################################################

# Check if random node is compatible with the rest of preselected nodes:
# 1. No ancestor
# 2. No descendant
# 3. No sister or sister descendant

is_compatible <- function(previous_nodes, node, tree) {
  
  # Check if previously selected nodes are ancestors of random node
  node_ancestors <- phangorn::Ancestors(tree, node, "all")
  if (any(previous_nodes %in% node_ancestors)) {
    return(FALSE) # Nodes are ancestors
  }
  
  # Check if previously selected nodes are descendants of random node
  node_descendants <- phangorn::Descendants(tree, node, "all")
  if (any(previous_nodes %in% node_descendants)) {
    return(FALSE)
  }
  
  # Check if previously selected nodes are siblings or sibling descendants
  node_siblings <- phangorn::Siblings(tree, node)
  if (any(previous_nodes %in% node_siblings)) {
    return(FALSE)
  }
  
  siblings_descendants <- unlist(phangorn::Descendants(tree, node_siblings, "all"))
  if (any(previous_nodes %in% siblings_descendants)) {
    return(FALSE)
  }
  
  return(TRUE)
}

select_nodes <- function(n_nodes, tree, selected_nodes) {
  
  while(length(selected_nodes) < n_nodes) {
    random_node <- sample(filtered_nodes$node, 1)
    
    if (is_compatible(selected_nodes, random_node, tree)) {
      selected_nodes <- c(selected_nodes, random_node)
    }
  }
  
  return(selected_nodes)
}

### MAIN ######################################################################

set.seed(777)

# Load data
tree <- ape::as.phylo(treeio::read.beast("../data/annotated_tree_noresis.nexus"))
load("../data/tree_tibble.rda")
load("../data/homoplasy_mutations.rda")
load("../data/n_node_tips.rda")
load("../data/n_sister_tips.rda")
load("../data/sisters.rda")

# Select mutations of interest (RoHO > 1)
top_RoHO_mutations <- homoplasy_nodes_annotated_byMutation %>%
  arrange(desc(RoHO)) %>%
  filter(RoHO > 1)

pvalues <- vector("list", nrow(top_RoHO_mutations))
n_nodes <- vector("list", nrow(top_RoHO_mutations))

for (top_mutation in 1:nrow(top_RoHO_mutations)) {
  
  homoplasy <- top_RoHO_mutations$mutation[[top_mutation]]
  
  load("../data/SNP_count.rda")
  
  n_appearances <- unname(snp_count[homoplasy]) 
  real_RoHO <- top_RoHO_mutations$RoHO[[top_mutation]]
  
  load("../data/homoplasy_nodes.rda")
  
  original_mutation_nodes <- homoplasy_nodes %>%
    filter(mutation == homoplasy) %>%
    left_join(tree_tibble, by = c("node", "label"))
  
  n_nodes[[top_mutation]] <- nrow(original_mutation_nodes)
  
  if (nrow(original_mutation_nodes) > 1) {
    root_distance_mean = mean(original_mutation_nodes$distance_to_root)
    root_distance_sd = sd(original_mutation_nodes$distance_to_root)
    
    top_limit <- root_distance_mean + 2 * root_distance_sd
    bottom_limit <- root_distance_mean - 2 * root_distance_sd
  } else {
    top_limit <- original_mutation_nodes$distance_to_root +  original_mutation_nodes$distance_to_root * 0.1
    bottom_limit <- original_mutation_nodes$distance_to_root - original_mutation_nodes$distance_to_root * 0.1
  }
  
 
  
  # Combine results back into tree_tibble
  filtered_nodes <- tree_tibble %>%
    mutate(
      n_node_tips = n_node_tips,
      sister = sisters,
      n_sister_tips = n_sister_tips
    ) %>%
    filter(n_node_tips >= 2,
           n_sister_tips >= 2,
           distance_to_root <= top_limit,
           distance_to_root >= bottom_limit)
  
  
  # Vector to store all the random generated RoHO values
  random_RoHO_values <- numeric(1000)
  
  for (i in 1:1000) {
    print(paste0("Random: ", i))
    
    # Vector to store selected nodes to add mutation
    # We start with a fresh preselected node each iteration
    selected_nodes <- c(sample(filtered_nodes$node, 1))
    
    mutation_nodes <- select_nodes(n_appearances, tree, selected_nodes)
    
    random_RoHO <- filtered_nodes %>%
      filter(node %in% mutation_nodes) %>%
      summarise(
        n_mut_alleles = sum(n_node_tips),
        n_wt_alleles = sum(n_sister_tips)
      ) %>%
      mutate(RoHO = n_mut_alleles / n_wt_alleles)
    
    random_RoHO_values[i] <- random_RoHO$RoHO
  }
  
  # Calculate mean and standard deviation
  mu <- mean(random_RoHO_values)
  sigma <- sd(random_RoHO_values)
  
  # Create the probability function using gaussfunc
  pvalue <- gaussfunc(real_RoHO, mu = mu, sigma = sigma)
  pvalues[[top_mutation]] <- pvalue
  
}

top_RoHO_mutations$pvalue <- unlist(pvalues)
top_RoHO_mutations$n_nodes <- unlist(n_nodes)

save(top_RoHO_mutations, file = "../data/top_RoHO_mutations.rda")
