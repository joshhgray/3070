"""
This code is adapted and modified from GitHub [Specifically, the mol_to_graph function]
Repo Owner: https://github.com/maxhodak
Script Author: https://github.com/dakoner
Fork: https://github.com/dakoner/keras-molecules
Script: https://github.com/dakoner/keras-molecules/blob/master/convert_rdkit_to_networkx.py

Licensed under MIT License for commerical use, modification, distribution, and private use.
(accessed on [6 February 2025])
"""

import networkx as nx
from rdkit import Chem, RDLogger

# Supress RDkit warnings
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
            mol_graph.add_node(atom.GetIdx(), 
                               element=atom.GetSymbol(),
                               atomic_num=atom.GetAtomicNum(),
                               formal_charge=atom.GetFormalCharge(),
                               chiral_tag=atom.GetChiralTag(),
                               hybridization=atom.GetHybridization(),
                               num_explicit_hs=atom.GetNumExplicitHs(),
                               is_aromatic=atom.GetIsAromatic(),
                               )
            
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