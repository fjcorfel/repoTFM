library(treeio)
library(ape)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)

# Read SNP table -> df
snp_table <- fread("../data/SNP_table_noresis.txt")

# Read splits -> vector
splits <- readLines("../data/mysplits_column.txt")
splits <- gsub("\\[|\\]", "", splits)
splits <- as.numeric(unlist(strsplit(splits, ", ")))

files <- list.files("../data/ancestral_results", pattern = "*.nexus",
                    full.names = TRUE, recursive = TRUE)

# Procesar cada uno de los archivos .nexus de ancestral_results/
result_tree <- NULL

# Fución para procesar las mutaciones de cada nodo
process_node_mutations <- function(node_mutations, n_file, snp_table, splits){
  if (is.null(node_mutations)) return(list(mut = NULL, ref = NULL))
  
  updated_mutations <- lapply(node_mutations, function(mutation){
    first_char <- str_extract(mutation, "^[A-Za-z-]")
    position <- as.numeric(str_extract(mutation, "\\d+"))
    last_char <- str_extract(mutation, "[A-Za-z-]$")
    
    # Si es el primer árbol, no se aumenta la posición respecto al split
    new_mutation <- mutation
    if (n_file != 1){
      position <- position + (splits[n_file])
      new_mutation <- (paste0(first_char, position, last_char)) 
    }
    
    # Obtener la posición de referencia de la mutación en la SNP_table
    ref_position <- snp_table$Position[position]
    
    # Se devuelve la mutación modificada y su posición de referencia
    list(new_mutation = new_mutation, ref_position = ref_position)
    
  })
  
  # Se devuelven todas las mutaciones modificadas y sus posiciones de referencia del nodo
  list(
    mut = sapply(updated_mutations, `[[`, "new_mutation"),
    ref = sapply(updated_mutations, `[[`, "ref_position")
    )
  }

for (n_file in seq_along(files)){
  
  tree <- treeio::read.beast(files[n_file])
  tree <- treeio::as_tibble(tree)
  
  # Crear plantilla inicial con el primer árbol
  if (n_file == 1) {
    result_tree <- tree
    result_tree$ref_mutation_position <- vector("list", nrow(result_tree))
    result_tree$n_mutations <- vector("list", nrow(result_tree))
  }
  
  # Procesar cada uno de los nodos
  processed_nodes <- lapply(seq_along(tree$mutations), function(n_node) {
    process_node_mutations(tree$mutations[[n_node]], n_file, snp_table, splits)
  })
  
  # Actualizar columnas en result_tree
  # Actualizar mutaciones
  result_tree$mutations <- lapply(seq_along(processed_nodes), function(n_node) {
    c(result_tree$mutations[[n_node]], processed_nodes[[n_node]]$mut)
  })
  
  # Actulizar posiciones de referencia de las mutaciones
  result_tree$ref_mutation_position <- lapply(seq_along(processed_nodes),
                                              function(n_node) {
    c(result_tree$ref_mutation_position[[n_node]],
      processed_nodes[[n_node]]$ref)
  })
}

# Eliminar duplicados de las mutaciones
result_tree$mutations <- lapply(result_tree$mutations, unique)

# Calcular el conteo de mutaciones después de eliminar duplicados
result_tree$n_mutations <- sapply(result_tree$mutations, length)

# Guardar el resultado final
save(result_tree, file = "parsed_ancestral_result.rda")

# Representar correlación entre longitud de ramas y número de mutaciones
ggplot(result_tree, aes(x = branch.length, y = n_mutations)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Correlación entre longitud de rama y número de mutaciones",
       x = "Longitud de rama",
       y = "Número de mutaciones") +
  theme_minimal()

# Guardar la gráfica
ggsave("correlation_plot.png", width = 10, height = 10)
