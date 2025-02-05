# Load required libraries
library(ape)
library(treeio)
library(parallel)  # For mclapply

# Read the tree and normalize branch lengths
tree <- treeio::read.newick("../data/tree.nwk")
tree$edge.length <- tree$edge.length / mean(tree$edge.length)

# Calculate parameters
T <- sum(tree$edge.length)  # Total tree length
n <- 2211                  # Number of inferred homoplasy events
n_traits <- 100         # Number of traits to simulate
rate <- n / T              # Rate for simulation

# Parallel simulation using mclapply
simulated_traits <- mclapply(1:n_traits, function(i) {
  # Optional: print the iteration number (note that output order may vary in parallel)
  print(i)
  
  # Simulate trait evolution along the tree
  rTraitDisc(
    phy = tree,
    model = "ER",
    k = 2,
    rate = rate,
    states = c("0", "1"),      # "0" = absence; "1" = presence of the mutation
    freq = c(0.5, 0.5),        # Neutral: equal equilibrium frequencies
    ancestor = FALSE  
  )
}, mc.cores = 14)

save(simulated_traits, file = "simulated_traits.rda")

# library(ape)
# 
# # Simulate a random tree with 10 tips
# set.seed(42)
# tree <- rtree(10)
# 
# # Simulate discrete trait evolution with two states (0 = no mutation, 1 = mutation)
# trait <- rTraitDisc(tree, model = "ER", k = 2, states = c(0, 1))
# 
# # Identify mutation events (where state changes from 0 to 1)
# mutation_nodes <- which(trait == 1)
# first_mutations <- setdiff(mutation_nodes, tree$edge[,2][trait[tree$edge[,1]] == 1])
# 
# # Create a modified version where only first mutations are kept
# modified_trait <- rep(0, length(trait))  # Default all to no mutation
# names(modified_trait) <- names(trait)
# modified_trait[first_mutations] <- 1  # Keep only first mutation appearances
# 
# # Print results
# print(trait)             # Original simulated traits
# print(modified_trait)    # Modified traits with mutations only at origin