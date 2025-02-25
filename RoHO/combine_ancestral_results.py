import pandas as pd

ancestral_result_noresis = pd.read_csv("../data/ancestral_result_noresis_fixed.csv")
ancestral_result_resis = pd.read_csv("../data/ancestral_result_resis_fixed.csv")
ancestral_result_noresis.rename(columns={"ref_mutation_position": "mutations"}, inplace=True)
ancestral_result_resis.rename(columns={"ref_mutation_position": "mutations"}, inplace=True)

ancestral_result_merged = ancestral_result_noresis.merge(ancestral_result_resis,
                                                         on="node",
                                                         suffixes=("_1", "_2"))

def combine_mutations(row):
    mut1 = row['mutations_1'] if pd.notna(row['mutations_1']) else ""
    mut2 = row['mutations_2'] if pd.notna(row['mutations_2']) else ""

    combined = set(mut1.split('|')) | set(mut2.split('|'))
    combined.discard("")  # Eliminar cadenas vacías (si alguna de las columnas estaba vacía)

    return '|'.join(sorted(combined)), len(combined)


ancestral_result_merged[['mutations', 'n_mutations']] = ancestral_result_merged.apply(combine_mutations, axis=1, result_type='expand')

ancestral_result = ancestral_result_merged[["node", "mutations", "n_mutations"]]

ancestral_result.to_csv("../data/ancestral_result.csv", index=False)
