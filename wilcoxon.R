library(dplyr)

load("../data/homoplasy_mutations.rda")

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
