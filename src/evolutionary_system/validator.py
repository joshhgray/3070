'''
Validation functions
'''
from rdkit import Chem, RDLogger #type: ignore

'''
Function to validate the authenticity of a given smiles string
'''
def validate_smiles(smiles_list):
    
    valid_smiles = []
    
    # Temporarily suppress RDLogger
    # So the Terminal isn't spammed with SMILES parse errors for every invalid mol
    RDLogger.DisableLog('rdApp.*')
    
    for smiles_str in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles_str)
            if mol:
                valid_smiles.append(smiles_str)
                
        except Exception as e:
            print("Error parsing {smiles_str}: {e}")
            
    # Re-enable RDLogger
    RDLogger.EnableLog('rdApp.*')
    
    return valid_smiles