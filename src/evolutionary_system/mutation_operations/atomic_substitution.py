from rdkit import Chem
import random

# Dictionary of possible atomic substitutions
VALID_SUBSTITUTIONS = {
    6: [7,8,16], # Carbon to either Nitrogen, Oxygen, or Sulfur
    7: [6,8], # Nitrogen to either Carbon or Oxygen
    8: [6,7,16], # Oxygen to either Carbon, Nitrogen or Sulfur
}

def atomic_substitution(rw_mol):
    """
    Performs and atomic substiution on an RDKit molecule.

    :param mol: RDKit mol object
    :returns: Mutated RDKit mol object if successful else the original, unedited mol
    """
    original = rw_mol

    if not isinstance(rw_mol, Chem.RWMol):
        rw_mol = Chem.RWMol(rw_mol)

    # List of all atoms in molecule 
    candidate_atoms = [atom.GetIdx() for atom in rw_mol.GetAtoms() if atom.GetAtomicNum() in VALID_SUBSTITUTIONS]

    if not candidate_atoms:
        return original

    # Choose a random atom from the original mol
    target_atom_idx = random.choice(candidate_atoms)
    target_atom = rw_mol.GetAtomWithIdx(target_atom_idx)
    original_atomic_num = target_atom.GetAtomicNum()

    # Select a valid replacement atom given selected atom from original mol
    replacement_atomic_num = random.choice(VALID_SUBSTITUTIONS[original_atomic_num])

    # Perform the substitution
    target_atom.SetAtomicNum(replacement_atomic_num)

    # Validate and sanitize
    mutated_mol = rw_mol.GetMol()
    try:
        Chem.SanitizeMol(mutated_mol)

    # Return original mol if substitution has failed
    except Exception as e:
        return original
    
    return mutated_mol