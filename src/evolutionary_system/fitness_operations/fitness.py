# Base fitness function
from rdkit import Chem
from rdkit.Chem import QED, AllChem, rdFingerprintGenerator
from rdkit.DataStructs import TanimotoSimilarity
import numpy as np
import networkx as nx
from src.evolutionary_system.utils.nx_graph_to_mol import nx_graph_to_mol

def calculate_qed(compound):
    """
    Calculate QED score of a given SMILES string

    :param structure: the SMILES-notation chemical notation of the molecule
    :returns: QED score of the given structure
    """
    try:
        mol = nx_graph_to_mol(compound)
        if mol: # Valid compound
            return QED.qed(mol)
        else:
            print(f"Invalid compound: {compound}.")
            return 0.0

    except Exception as e:
        print(f"Error calculating QED for {compound}: {e}")
        return 0.0


def calculate_diversity(population):
    """
    Calculate population diversity based on molecular graph representations of eah
    compound within every BGC
    
    :param population: NetworkX DiGraph representing and holidng the entire population
    :return: Diversity score
    """

    # Initialize fingerprint generator
    # set radius to 2 for faster processing - don't need the detail of 3, same with fingerprint size
    fingerprint_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)

    # Extract every compound from the population of BGCs
    compounds = [
        compound["mol_graph"] for node in population.nodes
        if "compounds" in population.nodes[node]
        for compound in population.nodes[node]["compounds"]
        if compound.get("mol_graph") is not None
    ]

    # Convert to RDKit Mol Object for calculation
    rdkit_mols = [mol for mol in compounds if isinstance(mol, Chem.Mol)]

    # Generate molecular fingerprint
    fingerprints = [fingerprint_generator.GetCountFingerprints(rdkit_mol) for rdkit_mol in rdkit_mols]

    # Calculate Tanimoto Similarity (Jaccard Index)
    similarity_scores = [
        TanimotoSimilarity(fingerprints[i], fingerprints[j])
        for i in range(len(fingerprints))
        for j in range(i+1, len(fingerprints))
    ]

    avg_similarity = float(np.mean(similarity_scores) if similarity_scores else 0)
    population_diversity = 1 - avg_similarity

    return population_diversity

