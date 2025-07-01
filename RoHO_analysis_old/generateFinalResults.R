library(dplyr)
library(data.table)

# Load results
load("../data/homoplasy_mutations_resis_annotated.rda")
load("../data/random_permutation_results.rda")
load("../data/ranking_results.rda")
load("../data/tt_wilcox_results.rda")
snp_count <- fread("../data/SNP_count.csv")

ranking_results <- ranking_results %>%
  select(Mutation, pvalue, adj_pvalue_Bonf, adj_pvalue_BH) %>%
  rename(mutation = Mutation,
         ranking_pvalue = pvalue,
         ranking_adj_pvalue_Bonf = adj_pvalue_Bonf,
         ranking_adj_pvalue_BH = adj_pvalue_BH)

tt_wilcox_results <- tt_wilcox_results %>%
  select(mutation, tt_pvalue, tt_adj_pvalue_Bonf, tt_adj_pvalue_BH,
         wilcox_pvalue, wilcox_adj_pvalue_Bonf, wilcox_adj_pvalue_BH)

random_permutation_results <- random_permutation_results %>%
  rename(rd_perm_pvalue = pvalue)

rd_perm_adj_pvalue_Bonf <- p.adjust(random_permutation_results$rd_perm_pvalue, method = "bonferroni")
rd_perm_adj_pvalue_BH <- p.adjust(random_permutation_results$rd_perm_pvalue, method = "BH")
random_permutation_results <- random_permutation_results %>%
  select(mutation, rd_perm_pvalue) %>%
  mutate(rd_perm_adj_pvalue_Bonf = rd_perm_adj_pvalue_Bonf,
         rd_perm_adj_pvalue_BH = rd_perm_adj_pvalue_BH)

  
homoplasy_mutations_final <- homoplasy_mutations_resis_annotated %>%
  left_join(random_permutation_results, by = "mutation") %>%
  left_join(ranking_results, by = "mutation") %>%
  left_join(tt_wilcox_results, by = "mutation") %>%
  left_join(snp_count, by = "mutation") %>%
  rename(n_snps = count,
         resis_drug = drug,
         resis_condifence = confidence) %>%
  select(mutation,
         synonym,
         Rv_number,
         n_mut_alleles,
         n_wt_alleles,
         RoHO,
         n_snps,
         resis_drug,
         resis_condifence,
         rd_perm_pvalue,
         rd_perm_adj_pvalue_Bonf,
         rd_perm_adj_pvalue_BH,
         ranking_pvalue,
         ranking_adj_pvalue_Bonf,
         ranking_adj_pvalue_BH,
         wilcox_pvalue,
         wilcox_adj_pvalue_Bonf,
         wilcox_adj_pvalue_BH,
         tt_pvalue,
         tt_adj_pvalue_Bonf,
         tt_adj_pvalue_BH)

save(homoplasy_mutations_final, file = "../data/homoplasy_mutations_final.rda")
fwrite(homoplasy_mutations_final, file = "../data/homoplasy_mutations_final.csv")
