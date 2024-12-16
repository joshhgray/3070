# Base fitness function
from rdkit import Chem # type: ignore
from rdkit.Chem import QED # type: ignore

def calculate_qed(smiles_str):
    """
    Calculate QED score of a given SMILES string
    """
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        if mol: # Valid compound
            return QED.qed(mol)
        else: # Invalid compound
            return 0

    except Exception as e:
        print(f"Error calculating QED for {smiles_str}: {e}")
        return 0