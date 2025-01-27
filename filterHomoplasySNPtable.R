library(data.table)
library(dplyr)
library(parallel)

# Load SNP table
snp_table <- fread("../data/SNP_table_noresis_mutations.txt")
# Vector of all SNPs
snps <- snp_table$Mutation
rm(snp_table)

# Load mutation table (ancestral_result -> result_tree)
load("../data/ancestral_result.rda")
result_tree <- result_tree %>%
  select(node, label, ref_mutation_position)
# List with mutations for each node
muts <- result_tree$ref_mutation_position
rm(result_tree)

# Create binary table (rows = nodes, cols = SNPs)
# 
binary_table <- mclapply(muts, function(node_muts) {

  # Check if SNP is in node mutations
  as.numeric(snps %in% node_muts)
}, mc.cores = 14)

binary_table <- as.data.frame(do.call(rbind, binary_table))

# Sum of the total of nodes containing each mutation
snp_count <- colSums(binary_table)
names(snp_count) <- names(snps)
save(snp_count, file = "../data/SNP_mutation_count.rda")

