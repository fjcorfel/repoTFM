from pathlib import Path
import pandas as pd
import numpy as np

# Main directory with population based results
main_dir = Path('../data/POPULATION_BASED')

# SNP table for annotation
snp_table = pd.read_csv('../data/SNP_table_final_redundant.txt',
                        sep='\t')

snp_table['Synonym'] = snp_table['Synonym'].replace(['', '-'], pd.NA)
snp_table = snp_table[['Position', 'Rv_number', 'Synonym']].rename(columns={
    'Synonym': 'synonym'
})

# Empty list for storing population based results
results = []

# Iterate over population directories
for population_dir in main_dir.iterdir():
    if population_dir.is_dir():
        csv_files = list(population_dir.glob('*agefilter40.csv'))
        
        if csv_files:
            csv_file = csv_files[0]
            
            population_results = pd.read_csv(csv_file)
            population_results['Population'] = population_dir.name
            results.append(population_results)
            
# Combine all the population based results into a single dataframe
combined_results = pd.concat(results, ignore_index=True)
combined_results = combined_results.drop(columns=['branch_length'])

# Annotate genes
combined_results['Position'] = combined_results['mutation'].str.extract(r'(\d+)').astype(int)
snp_table['Position'] = snp_table['Position'].astype(int)
snp_table = snp_table.drop_duplicates(subset='Position')

combined_results = pd.merge(
    combined_results,
    snp_table,
    on='Position',
    how='left'
)
combined_results = combined_results.drop(columns=['Position'])


# Sort results by adj pvalue
combined_results = combined_results.sort_values(by='ttest_adj_pvalue_BH', ascending=True)

# Save results
combined_results.to_csv('BL_results_agefilter40.csv', index=False)
            