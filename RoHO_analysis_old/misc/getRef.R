# Script for extracting ref sequence from SNP table

library(data.table)

snp_table <- fread("../data/SNP_table_noresis.txt")
header <- ">ref"
seq <- paste0(snp_table$WT, collapse = "")

ref_seq <- paste(header, seq, sep = "\n")

write(ref_seq, file = "ref.fas")
