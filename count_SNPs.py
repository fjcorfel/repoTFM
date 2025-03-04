import pandas as pd

# Load annotated tree (cleaned)
annotated_tree = pd.read_csv('../data/vietnam/annotated_tree_cleaned.csv')

# Extract mutations into a list
mutations = [m for sublist in annotated_tree['mutations'].str.split('|') for m in sublist]

# Count mutations 
mutation_counter = pd.Series(mutations).value_counts()

# Convert to DataFrame
mutation_counter_df = mutation_counter.reset_index()
mutation_counter_df.columns = ['mutation', 'count']

# Export results
mutation_counter_df.to_csv('../data/vietnam/SNP_count.csv', index=False)
