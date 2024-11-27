library(dplyr)

original_table <- read.table("SNP_table.txt", header = TRUE, sep = "\t")

filtered_table <- original_table %>%
  filter(is.na(Position_in_resistant_list))

write.table(filtered_table, file = "SNP_table_noresis.txt", sep = "\t", quote = FALSE, row.names = FALSE, na = "")
