import pandas as pd

DATASET = 'global'

# Load redundant SNP table for checking synonymous mutations
snp_table = pd.read_csv(f'../data/SNP_table_final_redundant.txt', sep='\t')
# Create complete mutation column
snp_table['mutation'] = snp_table['Position'].astype(str) + snp_table['ALT']
position_variants = snp_table[['mutation', 'Variant_type']]


# Load annotated tree (cleaned)
annotated_tree = pd.read_csv(f'../data/{DATASET}/annotated_tree_cleaned.csv')

# Extract mutations into a list
mutations = [m for sublist in annotated_tree['mutations'].str.split('|') for m in sublist]

# Count mutations 
mutation_counter = pd.Series(mutations).value_counts()

# Convert to DataFrame
mutation_counter_df = mutation_counter.reset_index()
mutation_counter_df.columns = ['mutation', 'count']

# Merge with position_variants to get variant type
mutation_counter_df = mutation_counter_df.merge(position_variants, on='mutation', how='left')

# Export results
mutation_counter_df.to_csv(f'../data/{DATASET}/SNP_count.csv', index=False)

# Count synonymous mutations in annotated tree
exploded_tree = annotated_tree.assign(mutation=annotated_tree['mutations'].str.split('|')).explode('mutation')
merged_exp_tree = exploded_tree.merge(mutation_counter_df, on='mutation', how='left')
synonym_count = merged_exp_tree[merged_exp_tree['Variant_type'] == 'synonymous_variant'].groupby('node')['mutation'].count().reset_index()
synonym_count.rename(columns={"mutation": "n_synonym_mutations"}, inplace=True)
annotated_tree = annotated_tree.merge(synonym_count, on='node', how='left').fillna(0)
annotated_tree['n_synonym_mutations'] = annotated_tree['n_synonym_mutations'].astype(int)

# Export annotated tree with synonymous mutation count
annotated_tree.to_csv(f'../data/{DATASET}/annotated_tree_cleaned.csv', index=False)
