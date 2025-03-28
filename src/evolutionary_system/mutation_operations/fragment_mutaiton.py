from rdkit import Chem
from rdkit.Chem import rdFMCS
import random


class FragmentMutation:
    def __init__(self, fragment_library):
        self.fragment_library = [Chem.MolFromSmiles(smiles) 
                                 for smiles in fragment_library
                                 if Chem.MolFromSmiles(smiles)]
        
    def apply(self, mol):
        if mol is None:
            return None

        # Pick a random fragment from the library
        # TODO - could make this probabilistic - based on fragment count numbers
        frag = random.choice(self.fragment_library)
        mutated_mol = Chem.CombineMols(mol, frag)

        # Failure = return original
        if Chem.MolToSmiles(mol) == Chem.MolToSmiles(mutated_mol):
            return mol # Return original
        
        # Confirm validity
        try:
            rw_mol = Chem.RWMol(mutated_mol)
            Chem.SanitizeMol(rw_mol)

            return rw_mol
        except Exception as e:
            print(f"Error when attempting fragment mutation: {e}")
            return None
