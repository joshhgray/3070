import random
from rdkit import Chem
from rdkit.Chem import rdFMCS, rdmolops
from src.evolutionary_system.crossover_operations.crossover import Crossover

class MCSCrossover(Crossover):
    """
    Uses Maximum Common Substructure (MSC) algorithm to identify largest common substructure 
    between two parent molecules, randomly chooses one and attaches a left over fragment from
    each of the parents.
    """
    def apply(self, parent1, parent2):
        # Sometimes overly time consuming - set timeout to 1 second
        # if nothing is found in 1 second, then MCS crossover has failed
        mcs = rdFMCS.FindMCS([parent1, parent2], timeout=2)
        mcs_smarts = Chem.MolFromSmarts(mcs.smartsString)
        mcs_bond = [bond.GetIdx() for bond in mcs_smarts.GetBonds()]

        if not mcs_bond:
            #print(f"NO BONDS, bank is closed.")
            return None

        # Break up each parent into frags by removing the MCS
        parent1_frags = Chem.FragmentOnBonds(parent1, mcs_bond, addDummies=True)
        parent2_frags = Chem.FragmentOnBonds(parent2, mcs_bond, addDummies=True)

        try:
            Chem.SanitizeMol(parent1_frags)
            Chem.SanitizeMol(parent2_frags)
        except Exception as e:
            #print(f"Error fragmenting molecules.")
            return None

        # Ensure validity
        parent1_frags = Chem.GetMolFrags(parent1_frags, asMols=True)
        parent2_frags = Chem.GetMolFrags(parent2_frags, asMols=True)

        # Need at least 2 fragments per parent to perform crossover operation
        if len(parent1_frags) < 2 or len(parent2_frags) < 2:
            return None

        # Randonly choose a remaining fragment from each parent
        frag1 = random.choice(parent1_frags)
        frag2 = random.choice(parent2_frags)

        offspring = Chem.RWMol(mcs_smarts)
        
        # Combine fragments
        offspring = Chem.CombineMols(mcs_smarts, frag1)
        offspring = Chem.CombineMols(offspring, frag2)

        if offspring:
            try:
                Chem.SanitizeMol(offspring)

                # Prevent mols with dummy frags from getting into population
                dummy_atom = Chem.MolFromSmarts("[*]")
                offspring = Chem.DeleteSubstructs(offspring, dummy_atom)
                Chem.SanitizeMol(offspring)

                return offspring 
            except:
                return None
        return None