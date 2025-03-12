from rdkit import Chem
from rdkit.Chem import QED, Crippen, rdMolDescriptors, Lipinski
from src.evolutionary_system.utils.nx_graph_to_mol import nx_graph_to_mol
from moses.metrics.SA_Score.sascorer import calculateScore

def calculate_qed(compound):
    """
    Calculate QED (Quantitative Estimate of Drug-likeness) score of a compound.

    :param compound: NetworkX Graph Object representing molecular compound.
    :returns: QED score of the given compound, 0 if error.
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

def calculate_sa_score(compound):
    """
    Calculates Synthetic Accessiblity (SA) score of given compound.

    :param compound: NetworkX Graph object representing molecular compound.
    :returns: SA Score else 10.0 (maximum SA score - least likely to be synthesizable).
    """
    try:
        mol = nx_graph_to_mol(compound)
        Chem.SanitizeMol(mol)
        return calculateScore(mol) if mol else 10.0
    except Exception as e:
        return 10.0
    
def calculate_logp(compound):
    """
    Calculates LogP (partial coefficient) of given compound. TODO - remove ?

    :param compound: NetworkX Graph object representing molecular compound.
    :returns: LogP score, 0 if error.
    """
    try:
        mol = nx_graph_to_mol(compound)
        Chem.SanitizeMol(mol)
        return Crippen.MolLogP(mol) if mol else 0.0
    except Exception:
        return 0
    
def calculate_molecular_weight(compound):
    """
    Calculates Molecular Weight of given compound.

    :param compound: NetworkX Graph object representing molecular compound.
    :returns: Molecular Weight, 0 if error.
    """
    try:
        mol = nx_graph_to_mol(compound)
        Chem.SanitizeMol(mol)
        return rdMolDescriptors.CalcExactMolWt(mol) if mol else 0.0
    except Exception:
        return 0
    
def lipinski_score(compound):
    """
    Calculates a Lipinski Rule of Five score using RDKit built-in Lipinski module.
    The score is the normalized sum of rules broken

    Ro5:
    Molecular Weight (mw) <= 500 Da
    LogP <= 5
    Hydrogen Bond Donors (hbd) <= 5
    Hydrogen Bond Acceptors (hba) <= 10
    Rotatable Bonds (rb) <= 10

    :param compound: NetwrokX Graph object representing molecular compound.
    :returns: A normalized score (0-1)
    """
    mol = nx_graph_to_mol(compound)

    if not mol:
        return 0.0
    
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rb = Lipinski.NumRotatableBonds(mol)

    broken_rules = sum([mw > 500, logp > 5, hbd > 5, hba > 10, rb > 10])
 
    ro5_score = max(0, (broken_rules / 5))

    return ro5_score




