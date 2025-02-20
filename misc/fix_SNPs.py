import pandas as pd
import re

# Fix SNP count
snp_count = pd.read_csv("SNP_count.csv")
mutations = snp_count["mutation"].to_list()
genotype_mutations = [mutation[1:] for mutation in mutations]
snp_count["mutation"] = genotype_mutations

# Group by mutation and sum the count value
snp_count_fixed = snp_count.groupby("mutation")["count"].sum().reset_index()

# Save the fixed SNP count
snp_count_fixed.to_csv("SNP_count_fixed.csv", index=False)

# ----------------------------------------------------------------------------

# Fix ancestral results
ancestral_result = pd.read_csv("ancestral_result_resis.csv")

# Remove mutations with gaps ("-")
ancestral_result["ref_mutation_position"] = ancestral_result["ref_mutation_position"].apply(
    lambda x: "|".join(m for m in x.split("|") if "-" not in m) if isinstance(x, str) else x
)

# Replace empty string values with NaN
ancestral_result.replace("", pd.NA, inplace=True)

# Remove first nucleotide of each mutation
ancestral_result["ref_mutation_position"] = ancestral_result["ref_mutation_position"].apply(
    lambda x: "|".join(re.sub(r"^[A-Z\-]", "", mutation) for mutation in x.split("|")) if isinstance(x, str) else x
)

# Add a new column 'n' which is the count of mutations in 'ref_mutation_position'
ancestral_result["n_mutations"] = ancestral_result["ref_mutation_position"].apply(
    lambda x: len(x.split("|")) if isinstance(x, str) else 0
)

# Export fixed ancestral results
ancestral_result = ancestral_result[["parent", "node", "branch.length",
                                     "label", "ref_mutation_position",
                                     "n_mutations"]]

ancestral_result.to_csv("ancestral_result_resis_fixed.csv", index=False)