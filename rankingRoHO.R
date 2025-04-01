library(dplyr)
library(tidyr)
library(data.table)

# Read files
RoHO_100 <- fread("../data/global/global_RoHO_agefilter100.csv")
RoHO_40 <- fread("../data/global/global_RoHO_agefilter40.csv")
RoHO_100_homoplasies <- fread("../data/global/global_RoHO_homoplasies_agefilter100.csv")
RoHO_40_homoplasies <- fread("../data/global/global_RoHO_homoplasies_agefilter40.csv")
nosyn_sites <- fread("../data/global/Syn_nosyn_sites.tsv") %>%
  select(Rvnumber, NoSyn_sites) %>%
  rename(Rv_number = Rvnumber) %>%
  filter(!is.na(NoSyn_sites))

# Main function
analyze_rank_wilcoxon <- function(data) {
  
  # Order mutations by desc RoHO and assign ranking
  data <- data %>%
    select(mutation, synonym, Rv_number, RoHO) %>%
    left_join(nosyn_sites, by="Rv_number", relationship = "many-to-one") %>%
    filter(!is.na(NoSyn_sites))
  
  max_NoSyn_sites <- max(data$NoSyn_sites)
  
  data <- data %>%
    mutate(scaled_weigth = round(max_NoSyn_sites / NoSyn_sites)) %>%
    arrange(desc(RoHO)) %>%
    mutate(ranking = rank(-RoHO, ties.method = "average"))
  
  # Calculate pvalue
  rank_mutation <- function(Rv_val, data) {
    # Get mutation info
    rank_gene <- data %>% filter(Rv_number == Rv_val)
    
    # Get other mutations info
    rank_others <- data %>% filter(Rv_number != Rv_val)
    
    # Ponderate info based on scaled weight
    rank_values_mutation <- rep(rank_gene$ranking, rank_gene$scaled_weigth)
    rank_values_others <- rep(rank_others$ranking, rank_others$scaled_weigth)
    
    wilcox.test(rank_values_mutation, rank_values_others, alternative = "less")$p.value
  }
  
  # Apply Wilcox to each mutation
  # data <- data %>%
  #   group_by(mutation) %>%
  #   mutate(pvalue = rank_mutation(mutation, data)) %>%
  #   ungroup()
  
  # Adjust pvalues
  data <- data %>%
    rowwise() %>%
    mutate(pvalue = rank_mutation(Rv_number, data)) %>%
    ungroup() %>%
    mutate(adj_pvalue_BH = p.adjust(pvalue, method = "BH"),
           adj_pvalue_bonferroni = p.adjust(pvalue, method = "bonferroni")) %>%
    distinct(Rv_number, synonym, .keep_all = TRUE)
  
  return(data)
}

print("Starting processing...")
RoHO_100_homoplasies_results <- analyze_rank_wilcoxon(RoHO_100_homoplasies)
print("25%")
RoHO_40_homoplasies_results <- analyze_rank_wilcoxon(RoHO_40_homoplasies)
print("50%")
RoHO_100_results <- analyze_rank_wilcoxon(RoHO_100)
print("75%")
RoHO_40_results <- analyze_rank_wilcoxon(RoHO_40)

save(RoHO_100_homoplasies_results, file = "../data/global/ranking_results/100years_homoplasies.rda")
save(RoHO_40_homoplasies_results, file = "../data/global/ranking_results/40years_homoplasies.rda")
save(RoHO_100_results, file = "../data/global/ranking_results/100years.rda")
save(RoHO_40_results, file = "../data/global/ranking_results/40years.rda")

print("Processing complete!")
