from rdkit import Chem
from rdkit.Chem import QED
from src.evolutionary_system.utils.nx_graph_to_mol import nx_graph_to_mol

def calculate_qed(compound):
    """
    Calculate QED score of a compound

    :param compound: NetworkX Graph Object representing molecular compound
    :returns: QED score of the given compound
    """
    try:
        mol = nx_graph_to_mol(compound)
        Chem.SanitizeMol(mol)
        if mol: # Valid compound
            return QED.qed(mol)
        else:
            #print(f"Invalid compound: {compound}.")
            return 0.0

    except Exception as e:
        #print(f"Error calculating QED for {compound}: {e}")
        return 0.0
