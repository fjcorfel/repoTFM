library(data.table)
library(dplyr)

# Load nodes with homoplasy
load("../data/homoplasy_nodes.rda")

# Read SNP table (complete version)
snp_table <- as_tibble(fread("../data/SNP_table_noresis.txt"))  
snp_table <- snp_table %>%
  select(Position, WT, ALT, Synonym, Rv_number) %>%
  mutate(Mutation = paste0(WT, Position, ALT)) %>%
  select(Mutation, Synonym, Rv_number)

# Add gene names and Rv number 
homoplasy_nodes <- homoplasy_nodes %>%
  mutate(
    synonym = snp_table$Synonym[match(mutation, snp_table$Mutation)],
    Rv_number = snp_table$Rv_number[match(mutation, snp_table$Mutation)]
  ) %>%
  select(node,
         label,
         mutation,
         synonym,
         Rv_number,
         n_mut_alleles,
         n_wt_alleles,
         RoHO)

# Save results
save(homoplasy_nodes, file = "../data/homoplasy_nodes_genes.rda")
