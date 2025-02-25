library(dplyr)
load("../data/homoplasy_nodes.rda")
load("../data/tree_tibble.rda")

n_nodes <- vector("list", nrow(top_RoHO_mutations))

for (n_mutation in 1:nrow(homoplasy_nodes)) {
  print(n_mutation)
  
  homoplasy <- homoplasy_nodes$mutation[[n_mutation]]
  
  mutation_nodes <- homoplasy_nodes %>%
    filter(mutation == homoplasy) %>%
    left_join(tree_tibble, by = c("node", "label"))
  
  n_nodes[[n_mutation]] <- nrow(mutation_nodes)
}

homoplasy_nodes$n_nodes <- unlist(n_nodes)

n_homoplasy_nodes <- homoplasy_nodes %>%
  select(mutation, n_nodes) %>%
  distinct()

save(n_homoplasy_nodes, file = '../data/n_nodes_per_homoplasy.rda')
