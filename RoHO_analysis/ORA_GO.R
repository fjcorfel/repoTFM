library(data.table)
library(clusterProfiler)
library(tidyverse)
library(GO.db)
library(AnnotationDbi)
library(stringr)

# Load gff file
gff <- fread("../data/Mycobacterium_tuberculosis_H37Rv_gff_v5.gff") %>%
  as_tibble()

colnames(gff) <- c("seqname", "source", "feature", "start", "end",
                   "score", "strand", "frame", "attributes")

genes <- gff %>% filter(feature == "CDS")

# Extract genes info from gff
genes <- genes %>%
  mutate(
    Rv_number = str_extract(attributes, "Locus=[^;]+") %>% str_replace("Locus=", ""),
    synonym = str_extract(attributes, "Name=[^;]+") %>% str_replace("Name=", ""),
    go_term = str_extract(attributes, "Gene Ontology=[^;]+") %>% str_replace("Gene Ontology=", "")
    ) %>%
  dplyr::select(Rv_number, synonym, go_term) %>%
  drop_na(go_term) %>%
  separate_rows(go_term, sep = ",")
  
term2gene <- genes %>% dplyr::select(go_term, Rv_number)

term2name <- AnnotationDbi::select(GO.db, keys = term2gene$go_term, columns = "TERM", keytype = "GOID")

# ORA
files <- c("../data/global/global_RoHO_homoplasies_agefilter40_phoR_renamed.csv",
           "../data/global/global_RoHO_homoplasies_agefilter100_phoR_renamed.csv",
           "../data/global/global_RoHO_agefilter40_phoR_renamed.csv",
           "../data/global/global_RoHO_agefilter100_phoR_renamed.csv")

ORA_results <- list()

for (i in seq_along(files)) {
  data <- fread(files[i]) %>%
    arrange(desc(RoHO)) %>%
    head(500) %>%
    filter(RoHO > 2)
    
  
  gene_list <- unique(data$Rv_number)
  
  ORA_results[[i]] <- clusterProfiler::enricher(gene = gene_list,
                                           TERM2GENE = term2gene,
                                           TERM2NAME = term2name,
                                           pvalueCutoff = 0.05)
}


save(ORA_results, file = "../data/global/ORA_results.rda")
