library(dplyr)
library(tidyr)

# Ejecutar en parsing_results/
files <- list.files(path = ".", pattern = "\\.txt$", full.names = TRUE)

# Function to read table in correct format
read_table <- function(file){
  table <- read.delim(file, stringsAsFactors = FALSE)
  
  table$mutations <- strsplit(table$mutations, "/")
  table$ref_positions <- strsplit(table$ref_positions, "/")
  return(table)
}

all_tables <- lapply(files, read_table)

# Combine all tables
combined_tables <- Reduce(function(x, y){
  dplyr::full_join(x, y, by = c("parent", "node", "branch.length", "label"),
                   suffix = c(".x", ".y"))
}, all_tables)

# Combine mutations and positions
combined_result <- combined_tables %>%
  rowwise() %>%
  mutate(
    mutations = paste(unique(unlist(c(mutations.x, mutations.y))),
                      collapse = "/"),
    ref_positions = paste(unique(unlist(c(ref_positions.x, ref_positions.y))),
                          collapse = "/")
  ) %>%
  select(parent, node, branch.length, label, mutations, ref_positions)

combined_result <- as_tibble(combined_result)

# Write result table
write.table(combined_result, "annotated_mutations.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)

# Store result in R object
combined_result$mutations <- strsplit(combined_result$mutations, "/")
combined_result$ref_positions <- strsplit(combined_result$ref_positions, "/")
combined_result <- combined_result %>%
  mutate(
    mutations = lapply(mutations, function(x) if (length(x) == 0) NULL else x),
    ref_positions = lapply(ref_positions, function(x) if (length(x) == 0) NULL else x),
    mutation_number = lapply(mutations, length)
  )

save(combined_result, file = "annotated_mutations.Rda")
