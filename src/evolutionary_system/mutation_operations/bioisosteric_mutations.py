from rdkit import Chem
from rdkit.Chem import FunctionalGroups, AllChem
import random
from src.evolutionary_system.mutation_operations.mutation import Mutation
import numpy as np

# TODO - expand and ensure cheminformatic relevancy
BIOISOSTERIC_REPLACEMENTS = {
    "hydroxyl": [ # Target --> Replacement
        ("[OH]", "[SH]"), # Hydroxyl --> Thiol
        ("[OH]", "[NH2]"), # Hydroxyl --> Amine
        ("[OH]", "[OCH3]"), # Hydroxyl --> Methoxy
    ],
    "carboxyl": [
        ("C(=O)[OH]", "C(=O)[SH]"), # Carboxyl --> Thiocarboxyl
        ("C(=O)[OH]", "C(=O)[C]"), # Carboxyl --> Ester
    ],
    "amine": [
        ("[NH2]", "[OH]"), # Amine --> Hydroxyl
        ("[NH2]", "[NH]C=O"), # Amine --> Amide
        ("[NH2]", "[SH]"), # Amine --> Thiol
    ],
    "fluorine": [
        ("[F]", "[Cl]"), # Fluorine --> Chlorine
        ("[F]", "[Br]"), # Fluorine --> Bromine
        ("[F]", "[CF3]"), # Fluorine --> Trifluoromethyl
    ],
}

class BioisostericMutation(Mutation):
    """
    Applies a bioisosteric transformation to an rdkit Mol.
    Can both perform and revert bioisosteric group transformations.
    Add/Revert ratio defines the likelyhood of either ocurring 80% add 20% revert.
    """
    def __init__(self, add_revert_ratio=0.8):
        self.add_revert_ratio = add_revert_ratio
        self.mode = "add"

    def apply(self, rw_mol):
        # Backup
        original = Chem.Mol(rw_mol)
        successful_mols = []

        # Select Add/Revert
        if random.random() > self.add_revert_ratio:
            self.mode = "revert"

        # Find any valid replacement sites from original molecule
        for replacements in BIOISOSTERIC_REPLACEMENTS.values():
            for original, replacement in replacements:
                
                if self.mode == "revert":
                    # If Reversion is toggled, reverse reaction order
                    original, replacement = replacement, original

                try:
                    # All possible reactions
                    rxn = AllChem.ReactionFromSmarts(f"{original}>>{replacement}")
                    # All possible products
                    products = rxn.RunReactants((original,)) # expects tuple

                    # Ensure validity of new mols
                    for product in products:
                        new_mol = product[0]
                        try:
                            Chem.SanitizeMol(new_mol)
                            successful_mols.append(new_mol)
                        except:
                            continue
                except:
                    continue
            """
            All possible outcomes are stored and a random one is chosen.
            This adds some overhead, but provides a greater probability of a reaction occuring.
            """
            if successful_mols:
                return random.choice(successful_mols)    
            # No possible transformations available
            else:
                return original