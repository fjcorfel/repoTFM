import pandas as pd

# Load SNP counts
snp_count_noresis = pd.read_csv("../data/SNP_count_noresis.csv")
snp_count_resis = pd.read_csv("../data/SNP_count_resis.csv")

# Merge snp counts
snp_count = pd.concat([snp_count_noresis, snp_count_resis]).groupby("mutation", as_index=False).sum()

# Export new count as csv
snp_count.to_csv("../data/SNP_count.csv", index=False)
