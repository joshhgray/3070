from rdkit import Chem
import csv

def load_initial_population(file_path="../../smiles_pool.csv"):
    
    initial_population = []
    
    try:
        with open(file_path, mode="r") as file:
            for line in file:
                #.strip() added because the csv is a little messy
                smiles_str = line.split(',')[0].strip()
                if smiles_str:
                    initial_population.append(smiles_str)
    except FileNotFoundError:
        print(f"Error - Can't find {file_path}. Unable to load population.")
        return []
    
    # Remove any invalid strings that may have gotten through
    initial_population = [smiles_str for smiles_str in initial_population
                         if Chem.MolFromSmiles(smiles_str) is not None]
    
    # Remove molecules with less than 5* molecules
    # *temporary number - TODO - find a good cutoff size
    initial_population = [smiles_str for smiles_str in initial_population
                          if Chem.MolFromSmiles(smiles_str).GetNumAtoms() > 5]
    
    
    return initial_population

