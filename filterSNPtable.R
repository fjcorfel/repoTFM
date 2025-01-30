library(data.table)
library(dplyr)

# Load SNP count
load("../data/SNP_count.rda")
# Filter SNPs that appear more than once
snp_count_filtered <- snp_count[snp_count > 5]

# Load SNP table and select mutations
snp_table_mutations <- data.table::fread("../data/SNP_table_noresis.txt") %>%
  select(Position, WT, ALT) %>%
  mutate(Mutation = paste0(WT, Position, ALT)) %>%
  select(Mutation)

filtered_mutations <- snp_table_mutations$Mutation[snp_table_mutations$Mutation %in% names(snp_count_filtered)]
save(filtered_mutations, file = "../data/filtered_mutations.rda")