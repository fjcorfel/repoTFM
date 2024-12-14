library(treeio)
library(ape)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)

# Read SNP table
snp_table <- fread("../data/SNP_table_noresis.txt")

# Read splits positions
splits <- readLines("../data/mysplits_column.txt")
splits <- gsub("\\[|\\]", "", splits)
splits <- as.numeric(unlist(strsplit(splits, ", ")))

# List files of annotated trees in nexus format
files <- list.files("../data/ancestral_results", pattern = "*.nexus",
                    full.names = TRUE, recursive = TRUE)

# Create template for the result tree
result_tree <- NULL

# Function to process individual node mutations of each tree
process_node_mutations <- function(node_mutations, n_file) {
  if (is.null(node_mutations) || all(node_mutations == "" )) {
    return(list(mut = NULL, ref = NULL))
  } 
  
  # List of updated mutations for the inidividual node
  # Each mutation is processed at once
  updated_mutations <- lapply(node_mutations, function(mutation) {

    first_char <- str_extract(mutation, "^[A-Za-z-]")
    position <- as.numeric(str_extract(mutation, "\\d+"))
    last_char <- str_extract(mutation, "[A-Za-z-]$")
    
   # Increase the mutation position according to splits
   # The first tree is not affected by split position
    new_mutation <- mutation
    if (n_file != 1){
      position <- position + (splits[n_file])
      new_mutation <- paste0(first_char, position, last_char)
    }
    
    # Obtain the reference position from the SNP table
    ref_position <- snp_table$Position[position]
    ref_position <- paste0(first_char, ref_position, last_char)
    
    # Get updated mutation positon and ref position for the individual mutation
    list(new_mutation = new_mutation, ref_position = ref_position)
  })
  
  # Return updated node mutations and reference positions
  list(
    mut = sapply(updated_mutations, `[[`, "new_mutation"),
    ref = sapply(updated_mutations, `[[`, "ref_position")
  )
}

# Process each annotated tree at once
for (n_file in seq_along(files)){
  
  print(paste0("Processing tree Nº", n_file, "..."))
  # Read tree as tibble
  tree <- treeio::read.beast(files[n_file])
  tree <- treeio::as_tibble(tree)
  
  # Update NULL template using first tree
  if (n_file == 1) {
    result_tree <- tree
    result_tree$ref_mutation_position <- vector("list", nrow(result_tree))
    result_tree$n_mutations <- vector("list", nrow(result_tree))
  }
  
  # List of processed nodes
  processed_nodes <- lapply(seq_along(tree$mutations), function(n_node) {
    process_node_mutations(tree$mutations[[n_node]], n_file)
  })
  
  # Update result_tree mutations
  result_tree$mutations <- lapply(seq_along(processed_nodes), function(n_node) {
    c(result_tree$mutations[[n_node]], processed_nodes[[n_node]]$mut)
  })
  
  # Update result_tree ref_mutation_position
  result_tree$ref_mutation_position <- lapply(seq_along(processed_nodes),
    function(n_node) {
  c(result_tree$ref_mutation_position[[n_node]], processed_nodes[[n_node]]$ref)
  })
}

# Remove mutation duplicates
result_tree$mutations <- lapply(result_tree$mutations, unique)

# Calculate mutation count for each node
result_tree$n_mutations <- sapply(result_tree$mutations, length)

# Save final result
save(result_tree, file = "ancestral_result.rda")

# Represent correlation between branch length and number of mutations
# Filter gaps (indels)
result_tree_no_gaps <- result_tree %>% 
  rowwise() %>%
  mutate(
    mutations = list(mutations[!grepl("-", mutations)]),
    n_mutations = length(mutations)
  ) %>% 
  ungroup()

# Plot the correlation
ggplot(result_tree_no_gaps, aes(x = branch.length, y = n_mutations)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Branch Length Vs Number of Mutations",
       x = "Branch Length",
       y = "Number of Mutations") +
  theme_minimal()

# Save the plot
ggsave("correlation_plot.png")
