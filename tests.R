library(dplyr)
library(data.table)

load("../data/homoplasy_mutations.rda")
load("../data/homoplasy_nodes.rda")
snp_table <- fread("../data/SNP_table.txt")

### OLD STRAT ###

ordered_df <- homoplasy_nodes_annotated_byMutation %>%
  select(mutation, synonym ,Rv_number, RoHO) %>%
  arrange(desc(RoHO))

mutations <- ordered_df$mutation
synonyms <- ordered_df$synonym
genes <- ordered_df$Rv_number
ranks <- 1:length(genes)

wilcox_results <- sapply(genes, function(gene) {
  target_rank <- ranks[genes == gene] 
  other_ranks <- ranks[genes  != gene]
  
  test <- wilcox.test(target_rank, other_ranks, alternative = "less")
  return(test$p.value)
})

wilcox_results_df <- data.frame(Mutation = mutations,
                                Synonym = synonyms,
                                Rv_number = genes,
                                pvalue = wilcox_results)

### NEW STRAT ###

tt_results <- homoplasy_nodes_annotated_byMutation %>%
  group_by(Rv_number) %>%
  summarise(
    n_mutations = n(),
    pvalue = if (n_mutations >= 2) {
      tryCatch({
        t.test(n_mut_alleles, n_wt_alleles, paired = TRUE)$p.value
      }, error = function(e) NA)
    } else NA
  ) %>%
  mutate(
    synonym = snp_table$Synonym[match(Rv_number, snp_table$Rv_number)]
  ) %>%
  select(synonym, Rv_number, n_mutations, pvalue)


wilcox_results <- homoplasy_nodes_annotated_byMutation %>%
  group_by(Rv_number) %>%
  summarise(
    n_mutations = n(),
    pvalue = if (n_mutations >= 2) {
      tryCatch({
        wilcox.test(n_mut_alleles, n_wt_alleles, paired = TRUE)$p.value
      }, warning = function(w) NA)
    } else NA
  ) %>%
  mutate(
    synonym = snp_table$Synonym[match(Rv_number, snp_table$Rv_number)]
  ) %>%
  select(synonym, Rv_number, n_mutations, pvalue)

tt_results <- homoplasy_nodes %>%
  group_by(mutation) %>%
  summarise(
    n_nodes = n(),
    pvalue = if (n_nodes >= 2) {
      tryCatch({
        t.test(n_mut_alleles, n_wt_alleles, paired = TRUE)$p.value
      }, error = function(e) NA)
    } else NA
   )

wilcox_results <- homoplasy_nodes %>%
  group_by(mutation) %>%
  summarise(
    n_nodes = n(),
    pvalue = if (n_nodes >= 2) {
      tryCatch({
        wilcox.test(n_mut_alleles, n_wt_alleles, paired = TRUE)$p.value
      }, warning = function(w) NA)
    } else NA
  )

tt_results$adj_pvalue <- p.adjust(tt_results$pvalue, method = "bonferroni")
wilcox_results$adj_pvalue <- p.adjust(wilcox_results$pvalue, method = "bonferroni")


### RoHO - 1 STRAT ### 
homoplasy_mutations <- homoplasy_nodes_annotated_byMutation %>%
  mutate(norm_RoHO = RoHO - 1)

tt_results <- homoplasy_mutations %>%
  group_by(Rv_number) %>%
  summarise(
    n_mutations = n(),
    pvalue = if (n_mutations >= 2) {
      tryCatch({
        t.test(norm_RoHO, mu = 0)$p.value
      }, error = function(e) NA)
    } else NA
  ) %>%
  mutate(
    synonym = snp_table$Synonym[match(Rv_number, snp_table$Rv_number)]
  ) %>%
  select(synonym, Rv_number, n_mutations, pvalue)

tt_results_all <- t.test(homoplasy_mutations$norm_RoHO, mu = 0)$p.value


