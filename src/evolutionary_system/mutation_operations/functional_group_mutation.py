import random
from rdkit import Chem
from src.evolutionary_system.mutation_operations.mutation import Mutation

FUNCTIONAL_GROUPS = {
    "hydroxylation": {
        "atoms": [(8, 1), (1, 1)], # OH
        # O-H single bond
        "bonds": [(0, 1, Chem.BondType.SINGLE)]
    },
    "methylation": {
        "atoms": [(6, 3)], # -CH3
        # Attatch with single bond
        "bonds": [(0, 1, Chem.BondType.SINGLE)]
    },
    "amination": {
        "atoms": [(7, 2), (1, 1)], # -NH2
        # attach with single bond
        "bonds": [(0, 1, Chem.BondType.SINGLE)]
    },
    "fluorination": {
        "atoms": [(9, 0)], # -F
        # No bond
        "bonds": []
    },
    "carboxylation": {
        "atoms": [(6, 0), (8, 1), (8, 1)], # COO
        # C=O C-O
        "bonds": [(0, 1, Chem.BondType.DOUBLE), (0, 2, Chem.BondType.SINGLE)]
    },
}

class FunctionalGroupMutation(Mutation):
    """
    Applies one of the functional groups to a valid atom.
    """
    def __init__(self, mutation_type):
        self.mutation_type = mutation_type
        self.group_data = FUNCTIONAL_GROUPS[mutation_type]

    def apply(self, rw_mol):
        # Backup - in case of mutation failure
        original = rw_mol

        # TODO - temp fix
        # Convert to RWMOl if necessary (sometimes GA may pass a Mol instead)
        if not isinstance(rw_mol, Chem.RWMol):
            rw_mol = Chem.RWMol(rw_mol)

        # Find valid attachment site
        attachment_atoms = [
            atom.GetIdx() for atom in rw_mol.GetAtoms()
            if atom.GetAtomicNum() in [6, 7, 8, 16] # C, N, O, S
            and not atom.GetIsAromatic()
            and atom.GetExplicitValence() < atom.GetTotalValence()
        ]

        if not attachment_atoms:
            return original
        
        # Pick a random valid Carbon atom to mutate
        target_atom = random.choice(attachment_atoms)

        # Add functional group
        new_atom_idxs = []
        for atomic_num, valence in self.group_data["atoms"]:
            new_atom = Chem.Atom(atomic_num)
            idx = rw_mol.AddAtom(new_atom)
            new_atom_idxs.append(idx)

        if not new_atom_idxs:
            return original
        
        rw_mol.AddBond(target_atom, new_atom_idxs[0], Chem.BondType.SINGLE)

        for start, end, bond_type in self.group_data["bonds"]:
            # Filter out if invalid bond idx
            if start >= len(new_atom_idxs) or end >= len(new_atom_idxs):
                return original
            rw_mol.AddBond(new_atom_idxs[start], new_atom_idxs[end], bond_type)

        try:
            # Sanitization requires Mol type
            mutated_mol = rw_mol.GetMol()
            Chem.SanitizeMol(mutated_mol)
        except Exception as e:
            return original
        
        return mutated_mol