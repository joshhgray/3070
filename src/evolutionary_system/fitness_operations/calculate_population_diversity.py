from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import TanimotoSimilarity
import numpy as np
import networkx as nx
from src.evolutionary_system.utils.nx_graph_to_mol import nx_graph_to_mol

def calculate_population_diversity(population, sample_size=50):
    """
    Calculate population-level diversity (across each compound) based on 
    molecular graph representations of each compound within every BGC
    
    :param population: NetworkX DiGraph representing and holidng the entire population
    :param sampel_size: size of sample to pull for estimating population diversity.
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
    # Convert from nx.Graph to to mol (non-re-writable required for fingerprint extraction)
    mols = [nx_graph_to_mol(mol, return_rwmol=False) for mol in compounds]

    # Generate molecular fingerprint
    fingerprints = [fingerprint_generator.GetCountFingerprint(mol) for mol in mols]
    # Calculate Tanimoto Similarity (Jaccard Index)
    similarity_scores = [
        TanimotoSimilarity(fingerprints[i], fingerprints[j])
        for i in range(len(fingerprints))
        for j in range(i+1, len(fingerprints))
    ]

    avg_similarity = float(np.mean(similarity_scores) if similarity_scores else 0)
    population_diversity = 1 - avg_similarity
    print(f"Pop diversity: {population_diversity}")

    return population_diversity