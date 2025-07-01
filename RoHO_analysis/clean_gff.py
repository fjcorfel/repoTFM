import pandas as pd

gff = pd.read_csv('../data/Mycobacterium_tuberculosis_H37Rv_gff_v5.gff',
                  sep='\t',
                  comment='#',
                  header=None)

attributes = gff.iloc[:, 8].str.split(';', expand=True)