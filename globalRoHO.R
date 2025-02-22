library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(ggplot2)

load("../data/tree_tibble.rda")
result_tree_noresis <- fread("../data/ancestral_result_noresis_fixed.csv")
result_tree_resis <- fread("../data/ancestral_result_resis_fixed.csv")
snp_table <- as_tibble(fread("../data/SNP_table_final_redundant.txt"))


# Only with nodes with 4 < n_tips (samples) < 100
# Calculate n of descendants
load("../data/n_node_tips.rda")
load("../data/n_sister_tips.rda")
global_RoHO <- tree_tibble %>%
  mutate(n_node_tips = n_node_tips,
         n_sister_tips = n_sister_tips) %>%
  filter(n_node_tips > 4 & n_node_tips < 100) %>%
  mutate(RoHO = n_node_tips / n_sister_tips) %>%
  mutate(mutations = ifelse(
    result_tree_noresis$ref_mutation_position[match(node, result_tree_noresis$node)] != "",
    result_tree_noresis$ref_mutation_position[match(node, result_tree_noresis$node)],
    result_tree_resis$ref_mutation_position[match(node, result_tree_resis$node)]
  )) %>%
  filter(mutations != "") %>%
  separate_rows(mutations, sep = "\\|") %>%
  rename(mutation = mutations) %>%
  group_by(mutation) %>%
  summarise(
    n_node_tips = sum(n_node_tips),
    n_sister_tips = sum(n_sister_tips),
    RoHO = n_node_tips / n_sister_tips
  ) %>%
  mutate(
    synonym = snp_table$Synonym[match(str_sub(mutation, 1, -2), snp_table$Position)],
    Rv_number = snp_table$Rv_number[match(str_sub(mutation, 1, -2), snp_table$Position)]
  ) %>%
  select(
    mutation, synonym, Rv_number, n_node_tips, n_sister_tips, RoHO
  )

save(global_RoHO, file = "../data/global_RoHO.rda")
fwrite(global_RoHO, file = "../data/global_RoHO.csv")

histogram <- global_RoHO %>%
  ggplot(aes(x = log2(RoHO))) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.4) +
  ggtitle("Distribution of RoHO values")

histogram
