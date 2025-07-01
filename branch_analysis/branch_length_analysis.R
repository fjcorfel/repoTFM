
# Load libraries
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(phangorn)
  library(parallel)
  library(pbmcapply)
  library(stringr)
  library(ape)
  library(treeio)
})

N_CORES <- 14

args <- commandArgs(trailingOnly = TRUE)

DATASET <- args[1]
TOTAL_TREE_YEARS <- 4000
YEARS_TO_ANALYZE <- 40



### FUNCTIONS -----------------------------------------------------------------

get_node_mutations <- function(node) {
  if (!node %in% annotated_tree$node) return(character(0))
  
  node_mutations <- unlist(annotated_tree[[which(annotated_tree$node == node), "mutations"]])
  return(node_mutations)
}


check_reversions <- function(node, mutation) {
  mutation_position <- as.numeric(str_extract(mutation, "\\d+"))
  
  descendants <- phangorn::Descendants(tree, node, "all")
  descendant_mutations <- unlist(lapply(descendants, get_node_mutations), use.names = FALSE)
  
  matches <- str_detect(descendant_mutations, paste0(mutation_position, "(?=\\D)"))
  return(any(matches))
}

check_mutations_in_wt <- function(node, mutation) {
  siblings <- phangorn::Siblings(tree, node)
  sibling_descendants <- unlist(phangorn::Descendants(tree, siblings, "all"))
  sibling_descendant_mutations <- unlist(lapply(sibling_descendants, get_node_mutations), use.names = FALSE)
  
  return(mutation %in% sibling_descendant_mutations)
}

calculate_n_tips <- function(node) {
  tips <- unlist(phangorn::Descendants(tree, node, "tips"))
  return(length(tips))
}

calculate_n_sibling_tips <- function(node) {
  siblings <- unlist(phangorn::Siblings(tree, node))
  sibling_tips <- unlist(phangorn::Descendants(tree, siblings, "tips"))
  return(length(sibling_tips))
}

get_node_branch_length <- function(node) {
  node_branch_length <- nodes_branch_length[[node, "branch.length"]]
}

test_branch_lengths <- function(node) {
  node_tips <- unlist(phangorn::Descendants(tree, node, "tips"))
  siblings <- unlist(phangorn::Siblings(tree, node))
  sibling_tips <- unlist(phangorn::Descendants(tree, siblings, "tips"))
  
  mutant_branch_lengths <- sapply(node_tips, get_node_branch_length)
  wt_branch_lengths <- sapply(sibling_tips, get_node_branch_length)
  
  # This block handles constant data error in t.test
  # This error happens when there is not branch length variation
  results <- tryCatch({
      t.test(mutant_branch_lengths, wt_branch_lengths, alternative = "less")

  }, error = function(e) {
    return(NA)
  })
  
  if (inherits(results, "htest")) {
    return(results$p.value)
  } else {
    return(NA)
  }
}

get_ancestors <- function(node) {
  ancestors <- phangorn::Ancestors(tree, node)
  return(ancestors)
}  

get_descendants <- function(node) {
  descendants <- phangorn::Descendants(tree, node, "all")
  return(descendants)
}  

get_n_synonym_mutations <- function(node) {
  mutation_idx <- which(annotated_tree$node == node)
  if (length(mutation_idx) > 0) {
    return(as.numeric(annotated_tree$n_synonym_mutations[mutation_idx]))
  } else return(0)
}

### INPUTS --------------------------------------------------------------------

# Load tree (nwk)
tree <- ape::as.phylo(treeio::read.beast.newick(paste0("../data/POPULATION_BASED/", DATASET, "/", "tree.nwk"))) 

# Load annotated_tree_cleaned (mutations per node and accumulated synonym mutations)
annotated_tree <- fread(paste0("../data/POPULATION_BASED/", DATASET, "/", "annotated_tree_cleaned.csv")) %>%
  as_tibble()


  
print("Calculating accumulated mutations root")
accumulated_mutations_root <- unlist(mclapply(annotated_tree$node, function(n) {
  node_mutations <- get_n_synonym_mutations(n)
  ancestors <- get_ancestors(n)
  ancestors_mutations <- sum(sapply(ancestors, function(ancestor) get_n_synonym_mutations(ancestor))) 
  
  total_mutations <- node_mutations + ancestors_mutations
  return(total_mutations)
}, mc.cores = N_CORES, mc.preschedule = FALSE))

print("Calculating accumulated mutations descendants")
accumulated_mutations_descendants <- unlist(mclapply(annotated_tree$node, function(n) {
  node_mutations <- get_n_synonym_mutations(n)
  descendants <- get_descendants(n)
  descendants_mutations <- sum(sapply(descendants, function(descendant) get_n_synonym_mutations(descendant))) 
  
  total_mutations <- node_mutations + descendants_mutations
  return(total_mutations)
}, mc.cores = N_CORES, mc.preschedule = FALSE))
  



annotated_tree$accumulated_mutations_root <- accumulated_mutations_root
annotated_tree$accumulated_mutations_descendants <- accumulated_mutations_descendants

distances_to_root <- node.depth.edgelength(tree)
distances_to_root <- data.frame(node = 1:length(distances_to_root), distance_to_root = distances_to_root)
annotated_tree <- annotated_tree %>%
  left_join(distances_to_root, by="node")


mutations_total_tree_years <- max(annotated_tree$accumulated_mutations_root)
mutations_years_to_analize <- round((YEARS_TO_ANALYZE * mutations_total_tree_years) / TOTAL_TREE_YEARS)

annotated_tree <- annotated_tree %>%
  filter(accumulated_mutations_descendants <= mutations_years_to_analize)
  

# Load annotated_tree (branch length per node)
nodes_branch_length <- fread(paste0("../data/POPULATION_BASED/", DATASET, "/", "annotated_tree.csv")) %>%
  select(node, branch.length) %>%
  as_tibble()

# Select mutations
snp_count <- fread(paste0("../data/POPULATION_BASED/", DATASET, "/", "SNP_count.csv")) #%>%
  #filter(count >= 5)
mutations <- snp_count$mutation
# Load SNP table for gene annotation
snp_table <- fread(paste0("../data/POPULATION_BASED/", DATASET, "/", "SNP_table_final.txt"))

### MAIN PROCESSING -----------------------------------------------------------

# Iterate over mutations
process_mutation <- function(n_mutation) {

  mutation <- mutations[n_mutation]
  
  # Filter nodes that have correct tip number size
  # Filter nodes that have the mutation
  # Check if mutation and WT branches are pure
  nodes <- annotated_tree %>%
    rowwise() %>%
    filter(mutation %in% unlist(mutations)) %>%
    ungroup()
  
  # Verify that nodes are not empty
  if (nrow(nodes) == 0) {
    return(tibble(node=integer(), mutation=character(),branch_length=character()))
  }
  
  nodes <- nodes %>%
    rowwise() %>%
    mutate(
      has_reversion = check_reversions(node, mutation),
      has_mutation_in_wt = check_mutations_in_wt(node, mutation)
    ) %>%
    ungroup() %>%
    filter(
      !has_reversion,
      !has_mutation_in_wt
    ) 
  
  # Verify that nodes are not empty
  if (nrow(nodes) == 0) {
    return(tibble(node=integer(), mutation=character(),branch_length=character()))
  }
  
  # Hacer aquí el cálculo
  # Sacar las longitudes de tips mutantes y tips wt
  # Hacer t.test entre ambos vectores -> guardar pvalor
  filtered_nodes <- nodes %>%
    mutate(mutation = mutation) %>%
    rowwise() %>%
    mutate(ttest_pvalue = test_branch_lengths(node)) %>%
    ungroup() %>%
    mutate(ttest_adj_pvalue_Bonf = p.adjust(ttest_pvalue, method = "bonferroni"),
           ttest_adj_pvalue_BH = p.adjust(ttest_pvalue, method = "BH")) %>%
    select(node, mutation, ttest_pvalue, ttest_adj_pvalue_Bonf, ttest_adj_pvalue_BH)
  
  return(filtered_nodes)
  
}

print("Starting branch lenght analysis...")

final_nodes <- pbmclapply(seq_along(mutations), function(n_mutation) {
  tryCatch({
    process_mutation(n_mutation)
  }, error = function(e) {
    write(paste("Error en mutación:", mutations[n_mutation], "->", conditionMessage(e)),
          file = "error_log.txt", append = TRUE)
    return(NULL)
  })
  
}, mc.cores = N_CORES, mc.preschedule = FALSE)

final_nodes <- do.call(rbind, final_nodes)

fwrite(final_nodes, file = paste0("../data/POPULATION_BASED/", DATASET, "/", "BL_final_nodes_", DATASET, "_agefilter", YEARS_TO_ANALYZE,".csv"))


# Group by mutation
final_mutations <- final_nodes %>%
  drop_na() %>%
  group_by(mutation) %>%
  summarise(adj_pvalues = list(ttest_adj_pvalue_BH),
            n = n(),
            significant_ratio = sum(ttest_adj_pvalue_BH <= 0.05) / n()) %>%
  ungroup()

final_mutations <- final_mutations %>%
  mutate(
    synonym = snp_table$Synonym[match(str_sub(mutation, 1, -2), snp_table$Position)],
    Rv_number = snp_table$Rv_number[match(str_sub(mutation, 1, -2), snp_table$Position)]
  ) %>%
  mutate(synonym = na_if(synonym, "")) %>%
  mutate(synonym = na_if(synonym, "-"))


fwrite(final_mutations, file = paste0("../data/POPULATION_BASED/", DATASET, "/", "BL_final_mutations_", DATASET, "_agefilter",YEARS_TO_ANALYZE,".csv"))


# Group by gene
final_genes <- final_nodes %>%
  drop_na() %>%
  mutate(
    synonym = snp_table$Synonym[match(str_sub(mutation, 1, -2), snp_table$Position)],
    Rv_number = snp_table$Rv_number[match(str_sub(mutation, 1, -2), snp_table$Position)]
  ) %>%
  group_by(Rv_number, synonym) %>%
  summarise(
    adj_pvalues = list(ttest_adj_pvalue_BH),
    n = n(),
    significant_ratio = sum(ttest_adj_pvalue_BH <= 0.05) / n(),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate(synonym = na_if(synonym, "")) %>%
  mutate(synonym = na_if(synonym, "-"))


fwrite(final_genes, file = paste0("../data/POPULATION_BASED/", DATASET, "/", "BL_final_genes_", DATASET, "_agefilter",YEARS_TO_ANALYZE,".csv"))
print(paste("Branch length analysis completed for dataset: ", DATASET))
