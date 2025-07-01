library(readxl)
library(dplyr)
library(tidyr)
library(clusterProfiler)

annotation <- read_excel("../data/newSangerFile_forAnnotation_withDeeperCategories.xlsx", 1) %>%
  dplyr::rename(
    Rv_number = Gene,
    synonym = Symbol,
    function_annot = `Functional annotation`
  )

term2gene <- annotation %>%
  dplyr::select(function_annot, Rv_number) %>%
  distinct()

term2name <- annotation %>%
  dplyr::select(function_annot) %>%
  distinct() %>%
  mutate(TERM = function_annot) %>%
  dplyr::rename(functional_category = function_annot)

files <- c("../data/global/global_RoHO_homoplasies_agefilter40_phoR_renamed.csv",
           "../data/global/global_RoHO_homoplasies_agefilter100_phoR_renamed.csv",
           "../data/global/global_RoHO_agefilter40_phoR_renamed.csv",
           "../data/global/global_RoHO_agefilter100_phoR_renamed.csv")

ORA_results <- list()

# Realizar el análisis ORA con las nuevas anotaciones
for (i in seq_along(files)) {
  data <- fread(files[i]) %>%
    arrange(desc(RoHO)) %>%
    head(500) %>%
    filter(RoHO > 2)
  
  gene_list <- unique(data$Rv_number)
  
  ORA_results[[i]] <- enricher(
    gene = gene_list,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pvalueCutoff = 0.05
  )
}

# Guardar resultados
save(ORA_results, file = "../data/global/ORA_results_new.rda")


library(clusterProfiler)
library(enrichplot)

for (i in seq_along(ORA_results)) {
  result <- ORA_results[[i]]
  
  # Verificar si el objeto es válido y contiene resultados con p ajustado significativo
  if (is.null(result) || nrow(result@result) == 0) {
    message(paste("⚠️  ORA", i, "está vacío o no tiene resultados."))
    next
  }
  
  # Filtrar los términos con p.adjust < 0.05
  sig_terms <- result@result %>% dplyr::filter(p.adjust < 0.05)
  
  if (nrow(sig_terms) == 0) {
    message(paste("⚠️  ORA", i, "no tiene términos significativamente enriquecidos (p.adjust < 0.05)."))
    next
  }
  
  # Mostrar el gráfico
  print(dotplot(result, showCategory = 20, title = paste("ORA Result", i)) +
          ggplot2::theme_minimal())
}

}


