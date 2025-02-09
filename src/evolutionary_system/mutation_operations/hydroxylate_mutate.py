import random
from rdkit import Chem

def hydroxylate_mutate(rw_mol):
    """
    Simulates hydroxylation mutation by adding an OH (-OH) group to a Carbon atom.
    Ensures that the target atom has available valency and is not aromatic.

    :param rw_mol: RDKit RWMol (modifiable molecule object).
    :returns: Mutated RDKit Mol object if successful, else original molecule.
    """
    # Backup - in case of mutation failure
    original = rw_mol

    # TODO - temp fix
    # Convert to RWMol if necessary (sometimes GA may pass an Mol instead)
    if not isinstance(rw_mol, Chem.RWMol):
        rw_mol = Chem.RWMol(rw_mol)

    carbon_atoms = []
    # Filter out non-aromatic, and full-valency atoms
    for atom in rw_mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetIsAromatic(): 
            continue
        if atom.GetExplicitValence() < 4: 
            carbon_atoms.append(atom.GetIdx())

    if not carbon_atoms:
        return original 

    # Pick a random valid Carbon atom to mutate
    target_atom = random.choice(carbon_atoms)

    # Add hydroxyl (-OH) group
    O = Chem.Atom(8) # Oxygen
    O_idx = rw_mol.AddAtom(O)
    rw_mol.AddBond(target_atom, O_idx, Chem.BondType.SINGLE)

    # Sanitization requires Mol type
    mutated_mol = rw_mol.GetMol()

    # Remove explicit hydrogens to prevent valency errors
    # This ocasionally fails, so it's wraped in a try/except
    try:
        mutated_mol = Chem.RemoveHs(mutated_mol)
    except Chem.AtomValenceException as e:
        return original

    # Sanitize molecule after modification
    try:
        Chem.SanitizeMol(mutated_mol)
    except Exception as e:
        print(f"Sanitization failed after mutation: {e}")
        return original

    return mutated_mol
