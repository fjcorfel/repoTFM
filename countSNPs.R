library(data.table)
library(dplyr)
library(parallel)

# Load SNP table
snp_table <- fread("../data/SNP_table_resis.txt")

snp_table_mutations <- snp_table %>%
  mutate(Mutation = paste0(WT, Position, ALT)) %>%
  select(Mutation)

snps <- snp_table_mutations$Mutation


# Load mutation table (ancestral_result -> result_tree)
load("../data/ancestral_result_resis.rda")
result_tree_resis <- result_tree_resis %>%
  select(node, label, ref_mutation_position)
# List with mutations for each node
muts <- result_tree_resis$ref_mutation_position

# Set up cluster
n_cores <- 14
cl <- makeCluster(n_cores)

# Share objects to the cluster
clusterExport(cl, varlist = c("snps"))

# Create binary table (rows = nodes, cols = SNPs)
binary_table <- parLapply(cl, muts, function(node_muts){
  snps %in% node_muts
})

# Stop the cluster
stopCluster(cl)

# Combine results into a single data frame
binary_table <- as.data.frame(do.call(rbind, binary_table))

# Sum of the total of nodes containing each mutation
snp_count <- colSums(binary_table)
names(snp_count) <- snps
save(snp_count, file = "../data/SNP_count_resis.rda")
