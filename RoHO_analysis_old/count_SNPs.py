import pandas as pd

# Read node-mutations table
data = pd.read_csv("../data/ancestral_result_noresis.csv")

# Select columns of interest
data = data[["node", "ref_mutation_position"]]

# Drop nodes with no mutations
data = data.dropna()

# Remove mutations with gaps ("-")
data["ref_mutation_position"] = data["ref_mutation_position"].apply(
    lambda x: "|".join(m for m in x.split("|") if "-" not in m)
)

# Drop rows with empty ref_mutation_position
data = data[data["ref_mutation_position"] != ""]

# Extract mutations into a list
mutations = [m for sublist in data["ref_mutation_position"].str.split("|") for m in sublist]

# Count mutations
mutation_counter = pd.Series(mutations).value_counts()

# Convert to DataFrame
mutation_counter_df = mutation_counter.reset_index()
mutation_counter_df.columns = ["mutation", "count"]

# Export to CSV
mutation_counter_df.to_csv("../data/SNP_count_noresis.csv", index=False)
