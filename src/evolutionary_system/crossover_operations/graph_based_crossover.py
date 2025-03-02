from src.evolutionary_system.utils.nx_graph_to_mol import nx_graph_to_mol
from src.data_pipeline.mol_to_graph import mol_to_graph
from rdkit import Chem
import networkx as nx
import random

def graph_based_crossover(parent1, parent2):
    """
    Performs graph-based crossover on two parent molecules. Random subgraphs from each parent
    are extracted and attempted to be merged together in a chemically valid structure.

    :params parent1, parent2: NetworkX Graph representing each of the parent molecules
    :returns: NetworkX Graph representing the offspring molecule. If crossover fails, no
              offspring is returned.
    """
    mol1 = nx_graph_to_mol(parent1, return_rwmol=True)
    mol2 = nx_graph_to_mol(parent2, return_rwmol=True)

    if mol1 is None or mol1.GetNumBonds() == 0 or mol2 is None or mol2.GetNumBonds() == 0:
        return None
    
    # Get bond indicies
    mol1_bonds = [bond.GetIdx() for bond in mol1.GetBonds()]
    mol2_bonds = [bond.GetIdx() for bond in mol2.GetBonds()]
    if not mol1_bonds or not mol2_bonds:
        return None
    
    # Choose random location on molecule to break
    # TODO - see if there is a more biologically accurate way to do this
    break_bond1 = random.choice(mol1_bonds)
    break_bond2 = random.choice(mol2_bonds)

    # Fix frequent out of bounds error
    if break_bond1 >= mol1.GetNumBonds() or break_bond2 >= mol2.GetNumBonds():
        return None

    # Extract fragments from breakpoint
    mol1_frag = Chem.FragmentOnBonds(mol1, [break_bond1])
    mol2_frag = Chem.FragmentOnBonds(mol2, [break_bond2])

    # Merge new fragments together
    offspring = Chem.CombineMols(mol1_frag, mol2_frag)

    # Convert new offspring to graph
    offspring_graph = mol_to_graph(Chem.MolToSmiles(offspring))

    # Only return offspring if crossover is successful
    return offspring_graph if offspring_graph else None