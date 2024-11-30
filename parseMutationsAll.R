library(treeio)
library(ape)
library(dplyr)
library(stringr)
library(data.table)

# Read SNP table -> df
#snp_table <- fread("./data/SNP_table_noresis.txt")

# Read splits -> vector
splits <- readLines("./data/mysplits.txt")
splits <- gsub("\\[|\\]", "", splits)
splits <- as.numeric(unlist(strsplit(splits, ", ")))

# Read all the annotated trees -> list[tibble]
tree_files <- list.files(path = ".",
                         pattern = "\\.nexus$",
                         full.names = TRUE)

annotated_trees <- lapply(tree_files, function(tree_file){
  tree <- treeio::read.beast(tree_file)
  treeio::as_tibble(tree)
})

# Function to modify each mutation

modify_mutation <- function(mutation, increment){
  
  # Extract mutation components using regex
  first_char <- str_extract(mutation, "^[A-Za-z-]")
  position <- as.numeric(str_extract(mutation, "\\d+"))
  last_char <- str_extract(mutation, "[A-Za-z-]$")
  
  # Increment the position number and reconstruct the mutation
  paste0(first_char, position + increment, last_char)
  
}

# Access split and respective tree
for (tree_idx in seq_along(annotated_trees)){
  current_tree <- annotated_trees[[tree_idx]]
  current_split <- splits[tree_idx]
  
  if (current_split == 0) {
    next
  }
  
  for (mut_idx in seq_along(current_tree$mutations)){
    mutations <- current_tree$mutations[[mut_idx]]
    
    if (!is.null(mutations)){
      new_mutations <- sapply(mutations,
                          modify_mutation,
                          increment = (current_split - 1))
      
      new_mutations <- str_trim(new_mutations)
      new_mutations <- unname(new_mutations)
      
      # Debugging
      #print(paste0(mutations, " -> ", new_mutations))
      
      current_tree$mutations[[mut_idx]] <- unname(new_mutations)
      annotated_trees[[tree_idx]] <- current_tree
      
    }
  }
}

