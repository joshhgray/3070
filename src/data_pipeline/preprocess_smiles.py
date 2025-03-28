import pandas as pd
from rdkit import Chem
from tqdm import tqdm
import os

def preprocess_smiles(smiles_csv, output_file):
    # progress bar since this file is > 1.5M compounds 
    tqdm.pandas()

    df = pd.read_csv("smiles_pool.csv", sep=";", on_bad_lines="skip")

    # Drop null data
    df = df[df["Smiles"].notna()]

    # Sanitize & filter out invalids
    def is_valid_mol(smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            Chem.SanitizeMol(mol)
            return True
        except:
            return False

    # Add valid extracted compounds to new df, then CSV
    df["valid_smiles"] = df["Smiles"].progress_apply(is_valid_mol)
    df = df[df["valid_smiles"]]

    df[["Smiles"]].to_csv(output_file, index=False)

    print(f"Successfully loaded {df.shape[0]:,} SMILES to population pool")

if __name__ == "__main__":
    cd = os.path.dirname(os.path.abspath(__file__))
    smiles_csv = os.path.join(cd, "../../smiles_pool.csv")
    output_file = os.path.join(cd, "../../preprocessed_smiles.csv")
    
    preprocess_smiles(smiles_csv, output_file)
