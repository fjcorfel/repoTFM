library(treeio)
library(ape)
library(dplyr)
library(stringr)
library(data.table)

# Habría que adaptar la ejecución para los 200 árboles...
# parallel --link Rscript parseMutations.R ::: $(ls data/annotated_trees) ::: {1..3}

# Read SNP table -> df
snp_table <- fread("../data/SNP_table_noresis.txt")

# Read splits -> vector
splits <- readLines("../data/mysplits.txt")
splits <- gsub("\\[|\\]", "", splits)
splits <- as.numeric(unlist(strsplit(splits, ", ")))

# Read annotated tree -> tibble
# Command Line
args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]
tree_number <- as.numeric(args[2])

# Manual input
# tree_file <- "../data/run_alignment_no_resis.002.nexus"
# tree_number <- 2

tree <- treeio::read.beast(tree_file)
tree <- treeio::as_tibble(tree)


# Function to modify each mutation
modify_mutation <- function(mutation, increment){
  
  # Extract mutation components using regex
  first_char <- str_extract(mutation, "^[A-Za-z-]")
  position <- as.numeric(str_extract(mutation, "\\d+"))
  last_char <- str_extract(mutation, "[A-Za-z-]$")
  
  # Increment the position number and reconstruct the mutation
  paste0(first_char, position + increment, last_char)
  
}

# First tree doesn't need modifications
if (tree_number != 1) {
  
  for (mut_idx in seq_along(tree$mutations)){
    mutations <- tree$mutations[[mut_idx]]
    
    if (!is.null(mutations)){
      new_mutations <- sapply(mutations,
                              modify_mutation,
                              increment = splits[tree_number] - 1)
      
      new_mutations <- str_trim(new_mutations)
      
      tree$mutations[[mut_idx]] <- unname(new_mutations)
      
      
    }
  }
}

# Buscar la posición de la mutación (tree$mutations)
# Buscar esa posición en la SNP table (número de fila)
# Colocar en la mutación la posición de la referencia -> Coger snp_table$position de la fila correspondiente


