import pandas as pd

# Load annotated tree
annotated_tree = pd.read_csv('../data/annotated_tree.csv')

# Drop nodes without any mutations
annotated_tree = annotated_tree.dropna(subset=['mutations'])

# Remove mutations with gaps ("-")
annotated_tree['mutations'] = annotated_tree['mutations'].apply(
    lambda x: '|'.join(m for m in x.split('|') if '-' not in m)
)

# Drop nodes with empty mutations after removing gaps
annotated_tree = annotated_tree[annotated_tree['mutations'] != '']

# Remove mutations first nucleotide 
annotated_tree['mutations'] = annotated_tree['mutations'].apply(
    lambda x: '|'.join([m[1:] for m in x.split('|')])
)

# Export results
annotated_tree.to_csv('../data/annotated_tree_cleaned.csv', index=False)
