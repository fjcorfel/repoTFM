library(data.table)
library(dplyr)
library(stringr)


# Load mutations (resis)
load("../data/ancestral_result_resis.rda")

# Load SNP table (resis) for checking positin equivalence
snp_table_resis <- fread("../data/SNP_table_resis.txt")

# Update mutations
update_mutations <- function(mutation_list, snp_table) {
  lapply(mutation_list, function(mutations) {
    
    if (is.null(mutations) || length(mutations) == 0) {
      return(NULL)
    }
    
    unname(sapply(mutations, function(mutation) {
      
      n_position <- as.numeric(str_extract(mutation, "\\d+"))
      real_position <- snp_table_resis$Position[n_position]
      str_replace(mutation, as.character(n_position), as.character(real_position))
      
    }))
  })
}

result_tree_resis$ref_mutation_position <- update_mutations(result_tree_resis$mutations, snp_table_resis)

result_tree_resis$n_mutations <- sapply(result_tree_resis$mutations, length)

save(result_tree_resis, file = "../data/ancestral_result_resis.rda")