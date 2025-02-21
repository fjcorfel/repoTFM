library(dplyr)
library(data.table)

load("../data/homoplasy_mutations.rda")
load("../data/homoplasy_nodes.rda")
snp_table <- fread("../data/SNP_table_final_redundant.txt")

### Mutation Ranking ###

ordered_df <- homoplasy_mutations %>%
  select(mutation, synonym ,Rv_number, RoHO) %>%
  arrange(desc(RoHO))

# Remove mutations without annotation
ordered_df <- ordered_df %>% filter(!is.na(Rv_number))

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

ranking_results <- data.frame(Mutation = mutations,
                                Synonym = synonyms,
                                Rv_number = genes,
                                pvalue = wilcox_results)

ranking_results$adj_pvalue_Bonf <- p.adjust(ranking_results$pvalue,
                                            method = "bonferroni")
ranking_results$adj_pvalue_BH <- p.adjust(ranking_results$pvalue,
                                          method = "BH")

save(ranking_results, file = "../data/ranking_results.rda")


### Paper Strategy ###

tt_results <- homoplasy_mutations %>%
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


wilcox_results <- homoplasy_mutations %>%
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
        wilcox.test(n_mut_alleles, n_wt_alleles, paired = TRUE)$p.value
    } else NA
  )

tt_results$adj_pvalue_Bonf <- p.adjust(tt_results$pvalue, method = "bonferroni")
tt_results$adj_pvalue_BH <- p.adjust(tt_results$pvalue, method = "BH")
tt_results <- tt_results %>%
  rename(tt_pvalue = pvalue,
         tt_adj_pvalue_Bonf = adj_pvalue_Bonf,
         tt_adj_pvalue_BH = adj_pvalue_BH)



wilcox_results$adj_pvalue_Bonf <- p.adjust(wilcox_results$pvalue, method = "bonferroni")
wilcox_results$adj_pvalue_BH <- p.adjust(wilcox_results$pvalue, method = "BH")
wilcox_results <- wilcox_results %>%
  rename(wilcox_pvalue = pvalue,
         wilcox_adj_pvalue_Bonf = adj_pvalue_Bonf,
         wilcox_adj_pvalue_BH = adj_pvalue_BH) %>%
  select(wilcox_pvalue, wilcox_adj_pvalue_Bonf, wilcox_adj_pvalue_BH)

tt_wilcox_results <- bind_cols(tt_results, wilcox_results)

save(tt_wilcox_results, file = "../data/tt_wilcox_results.rda")
