import pandas as pd
import argparse


def translate_mutation_position(mutations: str, snp_table) -> str:
    """
    Translates alignment mutation positions into whole genome positions.

    Args:
        mutations (str): Mutations separated by '|' character
        snp_table (DataFrame): SNP table with mutation annotation

    Returns:
        str: Mutations with fixed mutation position
    """
    translated_mutations = []
    for mutation in mutations.split('|'):
        # Extract mutation position number
        position = int(''.join(filter(str.isdigit, mutation)))
        real_position = snp_table.at[position - 1, 'Position']
        translated_mutations.append(f'{real_position}{mutation[-1]}')
        
    return '|'.join(translated_mutations)


def main(dataset):
    # Load annotated tree
    annotated_tree = pd.read_csv(f'../data/POPULATION_BASED/{dataset}/annotated_tree.csv')

    # Drop nodes without any mutations
    annotated_tree = annotated_tree.dropna(subset=['mutations'])

    # Remove mutations with gaps ("-")
    annotated_tree['mutations'] = annotated_tree['mutations'].apply(
        lambda x: '|'.join(m for m in x.split('|') if '-' not in m)
    )

    # Drop nodes with empty mutations after removing gaps
    annotated_tree = annotated_tree[annotated_tree['mutations'] != '']

    # This two steps aren't necessary for global dataset

    # Remove mutations first nucleotide 
    annotated_tree['mutations'] = annotated_tree['mutations'].apply(
        lambda x: '|'.join([m[1:] for m in x.split('|')])
    )

    # Fix mutation position -> translate alignment position into whole genome
    snp_table = pd.read_csv(f'../data/POPULATION_BASED/{dataset}/SNP_table_final.txt', sep='\t')
    annotated_tree['mutations'] = annotated_tree['mutations'].apply(
        lambda x: translate_mutation_position(x, snp_table)
    )

    # Export results
    annotated_tree.to_csv(f'../data/POPULATION_BASED/{dataset}/annotated_tree_cleaned.csv', index=False)
    print(f'Dataset {dataset}: cleaned and saved as annotated_tree_cleaned.csv')
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Clean annotated_tree.csv')
    parser.add_argument('--dataset', required=True, help='Dataset name (e.g., "Valencia", "Malawi")')
    args = parser.parse_args()
    main(args.dataset)
