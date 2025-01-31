import networkx as nx
from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.*')

def is_valid_smiles(smiles_str):
    """
    Verify whether a given string represent a valid smiles string
    
    :param smiles_str: smiles string to be tested
    :returns: True if valid, False if not
    """
    if not isinstance(smiles_str, str):
        return False

    try:
        mol = Chem.MolFromSmiles(smiles_str)

        if mol is None:
            return False
        if mol.GetNumAtoms() == 0:
            return False
        if mol.GetNumBonds() == 0 and mol.GetNumAtoms > 1:
            return False
        
        return True

    except Exception as e:
        print(f"Smiles validation error: {e}")
        return False
    

def mol_to_graph(smiles_str):
    """
    Uses RDkit to extract molecular form from smiles representation, 
    then converts to networkx directed graph

    :param smiles_str: smiles notation string representing a single, given BGC.
    :returns: networkX directed graph representing the atomic structure of the BGC.
    """
    if not is_valid_smiles(smiles_str):
        return None
    
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        mol_graph = nx.Graph()

        # Extract Atoms
        for atom in mol.GetAtoms():
            mol_graph.add_node(atom.GetIdx(), element=atom.GetSymbol())

        # Extract Bonds
        for bond in mol.GetBonds():
            mol_graph.add_edge(
                bond.GetBeginAtomIdx(),
                bond.GetEndAtomIdx(),
                bond_type=str(bond.GetBondType())
            )
            
        return mol_graph
    
    except Exception as e:
        print(f"Error in mol_to_graph: {e}.")
        return None