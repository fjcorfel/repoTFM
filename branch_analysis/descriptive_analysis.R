library(dplyr)
library(tidyr)
library(data.table)

BL_results_agefilter40 <- fread('BL_results_agefilter40.csv') %>%
  select(node, mutation, ttest_adj_pvalue_BH, Population, Rv_number, synonym) %>%
  filter(!is.na(ttest_adj_pvalue_BH))

BL_results_agefilter100<- fread('BL_results_agefilter100.csv') %>%
  select(node, mutation, ttest_adj_pvalue_BH, Population, Rv_number, synonym) %>%
  filter(!is.na(ttest_adj_pvalue_BH))

summary_mutations_agefilter40 <- BL_results_agefilter40 %>%
  group_by(mutation) %>%
  summarise(
    n_datasets = n_distinct(Population),
    datasets = paste(unique(Population), collapse = ", "),
    n_datasets_significative = sum(ttest_adj_pvalue_BH < 0.05),
    datasets_significative =  paste(unique(Population[ttest_adj_pvalue_BH < 0.05]), collapse = ", ")
  ) %>%
  arrange(desc(n_datasets_significative))


summary_genes_agefilter40 <- BL_results_agefilter40 %>%
  group_by(Rv_number) %>%
  summarise(
    synonym = first(na.omit(synonym)),
    n_datasets = n_distinct(Population),
    datasets = paste(unique(Population), collapse = ", "),
    n_datasets_significative = sum(ttest_adj_pvalue_BH < 0.05),
    datasets_significative =  paste(unique(Population[ttest_adj_pvalue_BH < 0.05]), collapse = ", ")
  ) %>%
  arrange(desc(n_datasets_significative))


summary_mutations_agefilter100 <- BL_results_agefilter100 %>%
  group_by(mutation) %>%
  summarise(
    n_datasets = n_distinct(Population),
    datasets = paste(unique(Population), collapse = ", "),
    n_datasets_significative = sum(ttest_adj_pvalue_BH < 0.05),
    datasets_significative =  paste(unique(Population[ttest_adj_pvalue_BH < 0.05]), collapse = ", ")
  ) %>%
  arrange(desc(n_datasets_significative))


summary_genes_agefilter100 <- BL_results_agefilter100 %>%
  group_by(Rv_number) %>%
  summarise(
    synonym = first(na.omit(synonym)),
    n_datasets = n_distinct(Population),
    datasets = paste(unique(Population), collapse = ", "),
    n_datasets_significative = sum(ttest_adj_pvalue_BH < 0.05),
    datasets_significative =  paste(unique(Population[ttest_adj_pvalue_BH < 0.05]), collapse = ", ")
  ) %>%
  arrange(desc(n_datasets_significative))
