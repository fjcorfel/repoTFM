library(dplyr)
library(pbmcapply)
library(data.table)


N_TOP_GENES <- 100
N_PERMUTATIONS <- 1000

# Load data dinamically
load_rda_file <- function(file_path) {
  load(file_path)
  return(global_RoHO)
}

# List of .rda files
rda_files <- c("../data/global/global_RoHO_homoplasies_agefilter40.rda")

for (rda_file in rda_files) {
  message(paste("Processing file:", rda_file))
   
  global_RoHO <- load_rda_file(rda_file)
  
  global_RoHO <- global_RoHO %>% arrange(desc(RoHO))
  
  observed_top_genes <- global_RoHO[1:N_TOP_GENES, ]
  
  top_genes <- unique(observed_top_genes$Rv_number)
  
  
  process_gene <- function(gene) {
    message(paste0("Processing ", gene))
    
    observed_count <- sum(observed_top_genes$Rv_number == gene)

    if (sum(global_RoHO$Rv_number == gene) < 3) {
      return(list(pvalue = NA, observed_count = observed_count))
    }
    
    
    # Perform random permutations
    random_count <- numeric(N_PERMUTATIONS)
    
    for (i in 1:N_PERMUTATIONS) {
      permuted_df <- global_RoHO %>% 
        mutate(Rv_number = sample(Rv_number))
      
      permuted_top_genes <- permuted_df[1:N_TOP_GENES, ]
      random_count[i] <- sum(permuted_top_genes$Rv_number == gene)
    }
    
    
    hist(random_count)
    # Compute pvalue
    # mu debería ser entorno a 0 siempre...
    mu <- mean(random_count)
    sigma <- sd(random_count)
    pvalue <- pnorm(observed_count, mean = mu, sd = sigma, lower.tail = FALSE)

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
  
  output_file <- paste0("permutation_results",
                        tools::file_path_sans_ext(basename(rda_file)),
                        ".csv")
  
  fwrite(final_results, paste0("../data/global/permutation_results/", output_file))
}

