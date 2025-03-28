from rdkit import Chem
from rdkit.Chem import AllChem
import random

class DeletionMutation:
    def apply(self, mol):
        """
        Applies one of the possible deletion methods.
        """
        DELETION_METHODS = [
            self.delete_random_atom_or_bond,
            self.delete_leaf_fragment,
            self.delete_terminal_ring
        ]
        random.shuffle(DELETION_METHODS)

        for method in DELETION_METHODS:
            rw_mol = Chem.RWMol(mol)
            try:
                new_mol = method(rw_mol)
                try:
                    Chem.SanitizeMol(new_mol)
                    return new_mol
                except:
                    continue

            except Exception as e:
                print(f"Error applying Deletion Mutation: {e}")
                continue
        
        # Return original mol if deletion has failed.
        return mol
    
    def delete_random_atom_or_bond(self, rw_mol):
        """
        TODO
        """
        # Extract atoms
        atoms = list(rw_mol.GetAtoms())
        # Choose at random
        atom_to_remove = random.choice(atoms)
        try:
            rw_mol.RemoveAtom(atom_to_remove.getIdx())
            return rw_mol
        except:
            return rw_mol
        

