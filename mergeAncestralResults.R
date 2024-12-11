library(treeio)
library(ape)
library(dplyr)
library(stringr)
library(data.table)

# Read SNP table -> df
snp_table <- fread("./data/SNP_table_noresis.txt")

# Read splits -> vector
splits <- readLines("./data/mysplits_column.txt")
splits <- gsub("\\[|\\]", "", splits)
splits <- as.numeric(unlist(strsplit(splits, ", ")))

files <- list.files("./data/ancestral_results/", pattern = "*.nexus", full.names = TRUE)


# Iterar sobre los árboles *.nexus de ancestral_results/
for (n_file in 1:length(files)) {

  # Si es el primer árbol, usar como "plantilla" y copiar -> añadir $ref_position
  if (n_file == 1){
    # Crear plantilla
    result_tree <- treeio::read.beast(files[n_file])
    result_tree <- treeio::as_tibble(result_tree)
    
    # Añadir columna ref_position
    result_tree$ref_mutation_position <- vector("list", nrow(result_tree))
    # Añadir columna con número de mutaciones
    result_tree$n_mutations <- vector("list", nrow(result_tree))
  }
  
  # A partir del segundo árbol -> realizar modificaciones
  else {
    # Leer árbol
    tree <- treeio::read.beast(files[n_file])
    tree <- treeio::as_tibble(tree)
    
    # Extraer lista de mutaciones
    tree_mutations_list <- tree$mutations
    
    # Recorrer cada uno de los nodos
    for (n_node in 1:length(tree_mutations_list)){
      node_mutations <- tree_mutations_list[[n_node]]
      
      # Si la fila no es NULL -> se hacen modificaciones
      if (!is.null(node_mutations)){
        
        # Recorrer cada una de las mutaciones de cada nodo
        for (n_mutation in 1:length(node_mutations)){
          # Extraer elementos
          first_char <- str_extract(node_mutations[n_mutation], "^[A-Za-z-]")
          position <- as.numeric(str_extract(node_mutations[n_mutation], "\\d+"))
          last_char <- str_extract(node_mutations[n_mutation], "[A-Za-z-]$")
          
          # Aumentar la posición respect al split
          position <- position + (splits[n_file] - 1)
          new_mutation <- paste0(first_char, position, last_char)
          
          # Extraer la posición de referencia de la mutación de la SNP_table
          ref_position <- snp_table$Position[position]
        
          # Añadir vector de mutaciones y ref_positions a la plantilla (append)
          result_tree$mutations[[n_node]] <- c(result_tree$mutations[[n_node]],
                                             new_mutation)
          
          result_tree$ref_mutation_position[[n_node]] <- c(result_tree$ref_mutation_position[[n_node]],
                                                          ref_position)
          
        }
      }
      
      # Añadir número de mutaciones resultante a cada nodo
      result_tree$n_mutations[[n_node]] <- length(result_tree$mutations[[n_node]])
    }
  }
}

save(result_tree, file = "parsed_ancestral_result.rda")
