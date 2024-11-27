awk 'NR > 1 {print $3}' SNP_table.txt | tr -d '\n' > secuencia.fasta
echo ">ref" | cat - secuencia.fasta > ref.fasta
rm secuencia.fasta