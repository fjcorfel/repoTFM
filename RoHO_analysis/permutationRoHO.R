library(dplyr)
library(pbmcapply)
library(data.table)


N_TOP_GENES <- 100
N_PERMUTATIONS <- 1000

# SNP table for normalization of genes based on length
snp_table <- fread("../data/SNP_table_final_redundant.txt") %>%
  mutate(gene_length = Gene_end - Gene_start + 1) %>%
  select(Rv_number, gene_length) %>%
  distinct(Rv_number, .keep_all = TRUE)

# Load data dinamically
load_rda_file <- function(file_path) {
  global_RoHO <- fread(file_path) %>%
    mutate(synonym = na_if(synonym, ""))
  
  return(global_RoHO)
}

# List of .rda files
rda_files <- c("../data/global/global_RoHO_homoplasies_agefilter40_phoR_renamed.csv",
               "../data/global/global_RoHO_homoplasies_agefilter100_phoR_renamed.csv",
               "../data/global/global_RoHO_agefilter40_phoR_renamed.csv",
               "../data/global/global_RoHO_agefilter100_phoR_renamed.csv")

for (rda_file in rda_files) {
  message(paste("Processing file:", rda_file))
   
  global_RoHO <- load_rda_file(rda_file) %>%
    left_join(snp_table, by=c("Rv_number"), relationship = "many-to-one")
  
  global_RoHO <- global_RoHO %>% arrange(desc(RoHO)) 
  
  observed_top_genes <- global_RoHO[1:N_TOP_GENES, ]
  
  top_genes <- unique(observed_top_genes$Rv_number)
  
  
  process_gene <- function(gene) {
    message(paste0("Processing ", gene))
    
    if (gene == "Rv0758_EXT" || gene == "Rv0758_INT") {
      gene_length <- 1458
      
    } else {
      gene_length <- observed_top_genes %>% 
        filter(Rv_number == gene) %>% 
        pull(gene_length) %>% 
        first()
    }
    
    
    
    observed_count <- sum(observed_top_genes$Rv_number == gene)

    if (sum(global_RoHO$Rv_number == gene) < 2 || is.na(gene_length)) {
      return(list(pvalue = NA, observed_count = observed_count))
    }
    
    
    # Perform random permutations
    random_count <- numeric(N_PERMUTATIONS)
    
    for (i in 1:N_PERMUTATIONS) {
      permuted_df <- global_RoHO %>% 
        mutate(Rv_number = sample(Rv_number))
      
      permuted_top_genes <- permuted_df[1:N_TOP_GENES, ]
      random_count[i] <- sum(permuted_top_genes$Rv_number == gene) / gene_length
    }
    
    # Compute pvalue
    mu <- mean(random_count)
    sigma <- sd(random_count)
    pvalue <- pnorm(observed_count / gene_length,
                    mean = mu,
                    sd = sigma,
                    lower.tail = FALSE)

    return(list(pvalue = pvalue, observed_count = observed_count))
    
    
  }
  
  set.seed(777)
  results <- pbmclapply(top_genes,
                        process_gene,
                        mc.cores = 12,
                        mc.preschedule = FALSE,
                        mc.set.seed = TRUE)
  
  final_results <- data.frame(
    Rv_number = top_genes,
    observed_top_count = sapply(results, function(x) x$observed_count),
    pvalue = sapply(results, function(x) x$pvalue)
  ) %>%
    left_join(global_RoHO %>% select(Rv_number, synonym),
              by="Rv_number") %>%
    distinct() %>%
    arrange(pvalue) %>%
    mutate(adj_pvalue_BH = p.adjust(pvalue, method = "BH")) %>%
    select(Rv_number, synonym, pvalue, adj_pvalue_BH, observed_top_count)
  
  output_file <- paste0("permutation_results_",
                        tools::file_path_sans_ext(basename(rda_file)),
                        ".csv")
  
  fwrite(final_results, paste0("../data/global/permutation_results/", output_file))
}

