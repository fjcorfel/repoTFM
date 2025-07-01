library(dplyr)
library(metap)
library(data.table)

BL_results_agefilter40 <- fread('BL_results_agefilter40.csv') %>%
  filter(!is.na(ttest_adj_pvalue_BH))
BL_results_agefilter100 <- fread('BL_results_agefilter100.csv') %>%
  filter(!is.na(ttest_adj_pvalue_BH))


mutation_combined_agefilter40 <- BL_results_agefilter40 %>%
  group_by(mutation) %>%
  summarise(
    combined_p_mutation = if(n() > 1) metap::sumlog(ttest_adj_pvalue_BH)$p else NA_real_
  ) %>% 
  filter(!is.na(combined_p_mutation)) %>%
  arrange(combined_p_mutation)

gene_combined_agefilter40 <- BL_results_agefilter40 %>%
  group_by(Rv_number) %>%
  summarise(
    combined_p_gene = if(n() > 1) metap::sumlog(ttest_adj_pvalue_BH)$p else NA_real_
  ) %>%
  filter(!is.na(combined_p_gene)) %>%
  arrange(combined_p_gene)

BL_results_agefilter40 <- BL_results_agefilter40 %>%
  left_join(mutation_combined_agefilter40, by = "mutation")

BL_results_agefilter40 <- BL_results_agefilter40 %>%
  left_join(gene_combined_agefilter40, by = "Rv_number")


mutation_combined_agefilter100 <- BL_results_agefilter100 %>%
  group_by(mutation) %>%
  summarise(
    combined_p_mutation = if(n() > 1) metap::sumlog(ttest_adj_pvalue_BH)$p else NA_real_
  ) %>% 
  filter(!is.na(combined_p_mutation)) %>%
  arrange(combined_p_mutation)

gene_combined_agefilter100 <- BL_results_agefilter100 %>%
  group_by(Rv_number) %>%
  summarise(
    combined_p_gene = if(n() > 1) metap::sumlog(ttest_adj_pvalue_BH)$p else NA_real_
  ) %>%
  filter(!is.na(combined_p_gene)) %>%
  arrange(combined_p_gene)

BL_results_agefilter100 <- BL_results_agefilter100 %>%
  left_join(mutation_combined_agefilter40, by = "mutation")

BL_results_agefilter100 <- BL_results_agefilter100 %>%
  left_join(gene_combined_agefilter40, by = "Rv_number")

fwrite(BL_results_agefilter40, file = "BL_results_agefilter40_pcombined.csv")
fwrite(BL_results_agefilter100, file = "BL_results_agefilter100_pcombined.csv")
