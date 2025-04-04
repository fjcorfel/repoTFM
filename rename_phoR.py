import pandas as pd

# Load data
snp_table = pd.read_csv('../data/SNP_table_final_redundant.txt', sep='\t')
snp_table['Mutation'] = snp_table['Position'].astype(str) + snp_table['ALT'].astype(str)

RoHO_data = pd.read_csv('../data/global/global_RoHO_agefilter100.csv')

ext_region = range(59, 177)

def update_gen(row):
    if row['synonym'] == 'phoR':
        mutation = row['mutation']
        aa_change = snp_table[snp_table['Mutation'] == mutation]['AA_change'].values
        
        if len(aa_change) > 0:
           aa_position = int(aa_change[0][1:-1])
           
           if aa_position in ext_region:
               new_synonym = 'phoR_EXT'
               new_Rv_number = 'Rv0758_EXT'
           else:
               new_synonym = 'phoR_INT'
               new_Rv_number = 'Rv0758_INT'
               
           return pd.Series([new_synonym, new_Rv_number])
           
    return pd.Series([row['synonym'], row['Rv_number']])

RoHO_data[['synonym', 'Rv_number']] = RoHO_data.apply(update_gen, axis=1)

RoHO_data.to_csv('../data/global/global_RoHO_agefilter100_phoR_renamed.csv',
                 index=False)

