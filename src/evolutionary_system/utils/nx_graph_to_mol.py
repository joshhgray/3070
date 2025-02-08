"""
This code is adapted and modified from GitHub 
Repo Owner: https://github.com/maxhodak
Script Author: https://github.com/dakoner
Fork: https://github.com/dakoner/keras-molecules
Script: https://github.com/dakoner/keras-molecules/blob/master/convert_rdkit_to_networkx.py

Licensed under MIT License for commerical use, modification, distribution, and private use.
(accessed on [6 February 2025])
"""
from rdkit import Chem
import networkx as nx

bond_type_map = {
    "SINGLE": Chem.BondType.SINGLE,
    "DOUBLE": Chem.BondType.DOUBLE,
    "TRIPLE": Chem.BondType.TRIPLE,
    "AROMATIC": Chem.BondType.AROMATIC
}

def nx_graph_to_mol(G, return_rwmol=True):
    """
    Converts a NetworkX molecular graph to an RDKit Mol object

    :param nx_graph: NetworkX Graph object representing molecule
    :returns: RDKit Mol object representing the same molecule (None if conversion fails) 
    """

    if not isinstance(G, nx.Graph):
        print(f"Error converting nx.Graph to rdkit mol for {G}")
        return None
    try:
        mol = Chem.RWMol()
        node_to_idx = {}

        atomic_nums = nx.get_node_attributes(G, 'atomic_num')
        #chiral_tags = nx.get_node_attributes(G, 'chiral_tag') TODO - see below
        formal_charges = nx.get_node_attributes(G, 'formal_charge')
        is_aromatic = nx.get_node_attributes(G, 'is_aromatic')
        hybridizations = nx.get_node_attributes(G, 'hybridization')
        num_explicit_hss = nx.get_node_attributes(G, 'num_explicit_hs')

        # Add atoms
        for node in G.nodes():
            try:
                a=Chem.Atom(atomic_nums[node])
                #a.SetChiralTag(chiral_tags[node]) - maybe causing errors TODO
                a.SetFormalCharge(formal_charges[node])
                a.SetIsAromatic(is_aromatic[node])
                a.SetHybridization(hybridizations[node])
                a.SetNumExplicitHs(num_explicit_hss[node])

                idx = mol.AddAtom(a)
                node_to_idx[node] = idx

            except Exception as e:
                print(f"Error processing atom {node}: {e}")

        # Add Bonds
        bond_types = nx.get_edge_attributes(G, 'bond_type')
        for edge in G.edges():
            try:
                first, second = edge
                if first not in node_to_idx or second not in node_to_idx:
                    print(f"Error: missing atoms for bonds {first}--{second}")
                    continue
                
                ifirst = node_to_idx[first]
                isecond = node_to_idx[second]

                # convert bond type to rdkit type
                bond_type = bond_types.get((first, second))
                if isinstance(bond_type, str):
                    bond_type = bond_type_map.get(bond_type)
                mol.AddBond(ifirst, isecond, bond_type)
            
            except Exception as e:
                print(f"Error adding bond {edge}: {e}")

        
        # Sanitize the mol - making sure it's chemically valid
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            print(f"Error sanitizing mol: {e}")
            return None

        return mol if return_rwmol else Chem.Mol(mol)

    except Exception as e:
        print(f"Error converting nx_graph to rdkit mol: {e}")
        return None