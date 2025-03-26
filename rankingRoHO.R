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
  rank_mutation <- function(mutation_val, data) {
    # Get mutation info
    rank_mutation <- data %>% filter(mutation == mutation_val)
    
    # Get other mutations info
    rank_others <- data %>% filter(mutation != mutation_val)
    
    # Ponderate info based on scaled weight
    rank_values_mutation <- rep(rank_mutation$ranking, rank_mutation$scaled_weigth)
    rank_values_others <- rep(rank_others$ranking, rank_others$scaled_weigth)
    
    wilcox.test(rank_values_mutation, rank_values_others, alternative = "less")$p.value
  }
  
  # Apply Wilcox to each mutation
  data <- data %>%
    group_by(mutation) %>%
    mutate(pvalue = rank_mutation(mutation, data)) %>%
    ungroup()
  
  # Adjust pvalues
  data <- data %>%
    mutate(adj_pvalue_BH = p.adjust(pvalue, method = "BH"),
           adj_pvalue_bonferroni = p.adjust(pvalue, method = "bonferroni"))
  
  return(data)
}

RoHO_100_homoplasies_results <- analyze_rank_wilcoxon(RoHO_100_homoplasies)
RoHO_40_homoplasies_results <- analyze_rank_wilcoxon(RoHO_40_homoplasies)
#RoHO_100_results <- analyze_rank_wilcoxon(RoHO_100)
#RoHO_40_results <- analyze_rank_wilcoxon(RoHO_40)

