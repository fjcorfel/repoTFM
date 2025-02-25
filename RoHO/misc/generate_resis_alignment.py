from Bio import SeqIO
import pandas as pd

# File name
alignment: str = "../data/alignment_resis_withref.fas"

# Load resistance idxs
# First idx is 1 -> list created on R
idxs = pd.read_csv("../data/resistance_idxs.csv")
idxs: list = idxs["idx"].values.tolist()

# Dict of samples {ID:seq}
samples:dict = {}

# Load file
for record in SeqIO.parse(alignment, "fasta"):
    samples[record.id] = record.seq
    
# Write new alignment
with open("alignment_resis_only.fas", "w") as file:
    for id, seq in samples.items():
        filtered_nucs:list = [seq[i - 1] for i in idxs]
        file.write(f">{id}\n{"".join(filtered_nucs)}\n") 
        