library(data.table)
library(dplyr)

snp_table <- data.table::fread("../data/SNP_table_noresis.txt") %>%
  select(Position, WT, ALT) %>%
  mutate(Mutation = paste0(WT, Position, ALT)) %>%
  select(Mutation)

write.table(snp_table, file = "../data/SNP_table_noresis_mutations.txt", quote = FALSE, row.names = FALSE)

tabla <- fread("../data/SNP_table_noresis_mutations.txt")
