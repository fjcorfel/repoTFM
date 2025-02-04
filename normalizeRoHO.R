library(dplyr)
library(data.table)

# Get gene lengths
snp_table <- fread("../data/SNP_table.txt")

gene_lengths <- snp_table %>%
  distinct(Rv_number, .keep_all = TRUE) %>%
  mutate(gene_length = Gene_end - Gene_start + 1) %>%
  select(Rv_number, gene_length)
  

# Calculate normalized RoHO value based on gene length
load("../data/homoplasy_mutations.rda")

homoplasy_nodes_annotated_byMutation_norm <- homoplasy_nodes_annotated_byMutation %>%
  mutate(
    gene_length = gene_lengths$gene_length[match(Rv_number, gene_lengths$Rv_number)]
  ) %>%
  mutate(
    norm_RoHO = RoHO / gene_length
  )

save(homoplasy_nodes_annotated_byMutation_norm, file = "../data/homoplasy_mutations_norm.rda")
