import pandas as pd
import subprocess
import io
import re
import glob
    

def annotate_mutation(mutation, res_files) -> dict:
    # Extract mutation position
    print(f"Annotating mutation {mutation}")
    mutation_position = re.findall(r"\d+", mutation)[0]

    # Grep cmd to find position across .res files
    grep_cmd = ["grep", "-wh", mutation_position] + res_files
    process = subprocess.run(grep_cmd, capture_output=True, text=True)
    
    # Convert output to df
    output = process.stdout
    if not output:
        return {"drug": "", "confidence": ""}
    
    df = pd.read_csv(io.StringIO(output), header=None).drop_duplicates(subset=[6, 7],
                                                                       keep="first")
    # Select columns of interest (Drug, Confidence)
    if df.shape[0] == 1:
        drug = str(df.iloc[:, 6].values[0])
        confidence = str(df.iloc[:, 7].values[0])
        
    else:
        drug = "; ".join(map(str, df.iloc[:, 6].to_list()))
        confidence = "; ".join(map(str, df.iloc[:, 7].to_list()))
        
    return {"drug": drug, "confidence": confidence}
    
    

def main() -> None:
    res_files: list = glob.glob("../data/res_files/*.res")
    if not res_files: 
        print("No .res found!")
    
    mutation_table = pd.read_csv("../data/homoplasy_mutations.csv")
    
    # Empty dict for storing mutation annotations
    drug_annotations = {}
    confidence_annotations = {}
    
    # Annotate each of the mutations
    for mutation in mutation_table["mutation"].to_list():
        annotation = annotate_mutation(mutation, res_files)
        
        drug_annotations[mutation] = annotation["drug"]
        confidence_annotations[mutation] = annotation["confidence"]
    
    mutation_table["drug"] = mutation_table["mutation"].map(drug_annotations)
    mutation_table["confidence"] = mutation_table["mutation"].map(confidence_annotations)
    
    mutation_table.to_csv("../data/homoplasy_mutations_resis_annotated.csv", index=False)
    
    print("\nProcess finished!")
    
    
if __name__ == "__main__":
    main()
