library(dplyr)
library(data.table)

# Load homoplasy nodes (resis and noresis)
load("../data/homoplasy_nodes_resis.rda")
load("../data/homoplasy_nodes_noresis.rda")

# Merge resis and noresis homoplasy nodes
homoplasy_nodes <- rbind(homoplasy_nodes, homoplasy_nodes_resis)
save(homoplasy_nodes, file = "../data/homoplasy_nodes.rda")

# Anotate mutations
print("Annotating mutations...")

# Read SNP table (complete version)
snp_table <- as_tibble(fread("../data/SNP_table.txt"))
snp_table <- snp_table %>%
  select(Position, WT, ALT, Synonym, Rv_number) %>%
  mutate(Mutation = paste0(WT, Position, ALT)) %>%
  select(Mutation, Synonym, Rv_number)

# Add gene names and Rv number 
homoplasy_nodes_annotated <- homoplasy_nodes %>%
  mutate(
    synonym = snp_table$Synonym[match(mutation, snp_table$Mutation)],
    Rv_number = snp_table$Rv_number[match(mutation, snp_table$Mutation)]
  ) %>%
  select(
    node,
    label,
    mutation,
    synonym,
    Rv_number,
    n_mut_alleles,
    n_wt_alleles,
    RoHO)

# Save results
save(homoplasy_nodes_annotated, file = "../data/homoplasy_nodes_annotated.rda")

# Group by mutation
homoplasy_nodes_annotated_byMutation <- homoplasy_nodes_annotated %>%
  group_by(mutation) %>%
  summarise(
    n_mut_alleles = sum(n_mut_alleles),
    n_wt_alleles = sum(n_wt_alleles),
    RoHO = n_mut_alleles / n_wt_alleles
  ) %>%
  mutate(
    synonym = snp_table$Synonym[match(mutation, snp_table$Mutation)],
    Rv_number = snp_table$Rv_number[match(mutation, snp_table$Mutation)]
  ) %>%
  select(
    mutation,
    synonym,
    Rv_number,
    n_mut_alleles,
    n_wt_alleles,
    RoHO
  ) %>%
  ungroup()

save(homoplasy_nodes_annotated_byMutation, file = "../data/homoplasy_mutations.rda")

# Group by gene
homoplasy_nodes_annotated_byGene <- homoplasy_nodes_annotated %>%
  group_by(Rv_number) %>%
  summarise(
    n_mut_alleles = sum(n_mut_alleles),
    n_wt_alleles = sum(n_wt_alleles),
    RoHO = n_mut_alleles / n_wt_alleles
  ) %>%
  mutate(
    synonym = snp_table$Synonym[match(Rv_number, snp_table$Rv_number)]
  ) %>%
  select(
    Rv_number,
    synonym,
    n_mut_alleles,
    n_wt_alleles,
    RoHO
  ) %>%
  ungroup()

save(homoplasy_nodes_annotated_byGene, file = "../data/homoplasy_genes.rda")
