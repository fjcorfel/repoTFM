suppressPackageStartupMessages({
  library(treeio)
  library(dplyr)
  library(data.table)
})

# Main directory with population based datasets
BASE_DIR="../data/POPULATION_BASED"

# Find all annotated_tree.nexus
nexus_files <- list.files(BASE_DIR,
                          pattern = "annotated_tree.nexus",
                          recursive = TRUE,
                          full.names = TRUE)

# Convert nexus to csv
convert_nexus_to_csv <- function(nexus_file) {
  
  # Read nexus file
  tree <- as_tibble(treeio::read.beast(nexus_file)) %>%
    mutate(mutations = lapply(mutations, function(x) if (is.null(x)) NA else x))
  
  # Create output name
  output_file <- gsub("treetime_results/annotated_tree.nexus",
                      "annotated_tree.csv",
                      nexus_file)
  
  # Write output csv
  fwrite(tree,
         file = output_file)
  
  message("Converted file: ", output_file)
  
  return(output_file)
}

output_files <- lapply(nexus_files, convert_nexus_to_csv)

message("All TreeTime '.nexus' files have been converted to '.csv'")
