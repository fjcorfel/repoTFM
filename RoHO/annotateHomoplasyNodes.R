library(dplyr)
library(data.table)
library(stringr)

# Load homoplasy nodes (resis and noresis)
load("../data/homoplasy_nodes_resis.rda")
load("../data/homoplasy_nodes_noresis.rda")

# Fix mutation column
homoplasy_nodes_noresis <- homoplasy_nodes_noresis %>%
  mutate(mutation = as.character(mutation$mutation))

homoplasy_nodes_resis <- homoplasy_nodes_resis %>%
  mutate(mutation = as.character(mutation$mutation))

# Merge resis and noresis homoplasy nodes
homoplasy_nodes <- bind_rows(homoplasy_nodes_noresis, homoplasy_nodes_resis)
save(homoplasy_nodes, file = "../data/homoplasy_nodes.rda")

# Anotate mutations
# Read SNP table (complete version)
snp_table <- as_tibble(fread("../data/SNP_table_final_redundant.txt"))
snp_table <- snp_table %>%
  select(Position, Synonym, Rv_number)

# Add gene names and Rv number 
homoplasy_nodes <- homoplasy_nodes %>%
  mutate(
    synonym = snp_table$Synonym[match(str_sub(mutation, 1, -2), snp_table$Position)],
    Rv_number = snp_table$Rv_number[match(str_sub(mutation, 1, -2), snp_table$Position)]
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

homoplasy_nodes <- homoplasy_nodes %>%
  mutate(synonym = na_if(synonym, "")) %>%
  mutate(synonym = na_if(synonym, "-"))

# Save results
save(homoplasy_nodes, file = "../data/homoplasy_nodes_annotated.rda")
fwrite(homoplasy_nodes, file = "../data/homoplasy_nodes_annotated.csv")

# Group by mutation
homoplasy_mutations <- homoplasy_nodes %>%
  group_by(mutation) %>%
  summarise(
    n_mut_alleles = sum(n_mut_alleles),
    n_wt_alleles = sum(n_wt_alleles),
    RoHO = n_mut_alleles / n_wt_alleles
  ) %>%
  mutate(
    synonym = snp_table$Synonym[match(str_sub(mutation, 1, -2), snp_table$Position)],
    Rv_number = snp_table$Rv_number[match(str_sub(mutation, 1, -2), snp_table$Position)]
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

homoplasy_mutations <- homoplasy_mutations %>%
  mutate(synonym = na_if(synonym, "")) %>%
  mutate(synonym = na_if(synonym, "-"))


save(homoplasy_mutations, file = "../data/homoplasy_mutations.rda")
fwrite(homoplasy_mutations, file = "../data/homoplasy_mutations.csv")


# Group by gene
homoplasy_genes <- homoplasy_nodes %>%
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

homoplasy_genes <- homoplasy_genes %>%
  mutate(synonym = na_if(synonym, "")) %>%
  mutate(synonym = na_if(synonym, "-"))

save(homoplasy_genes, file = "../data/homoplasy_genes.rda")
fwrite(homoplasy_genes, file = "../data/homoplasy_genes.csv")
