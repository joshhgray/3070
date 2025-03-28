from rdkit import DataStructs
from rdkit.Chem import rdMolDescriptors
import numpy as np
import random

def calculate_population_diversity(population, sample_size=50):
    """
    Calculate population-level diversity based on the Tanimoto Similarity 
    score (using the BulkTanimotoSimilarity function from RDKit) of a 20% 
    random sample of the population.
    
    :param population: NetworkX DiGraph representing and holidng the entire population
    :param sampel_size: size of sample to pull for estimating population diversity.
    :return: Diversity score
    """

    # Extract each of the active mols in the population
    mols = []
    for node, data in population.nodes(data=True):
        compounds = data.get("compounds")
        if compounds and compounds[0].get("mol") is not None:
            mols.append(compounds[0]["mol"])

    # Dynamically adjust sample size as 20% of current population size
    sample_size = int(0.2 * len(mols))

    if mols:
        sampled_mols = random.sample(mols, sample_size)
    else:
        return 0.0

    # Generate molecular fingerprint
    fingerprints = [rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024) for mol in sampled_mols]

    similarity_scores = []
    for i in range(len(fingerprints)):
        similarity_score = DataStructs.BulkTanimotoSimilarity(fingerprints[i], fingerprints[i+1:])
        similarity_scores.extend(similarity_score)

    if similarity_scores:
        avg_similarity = np.mean(similarity_scores)
    else:
        return 0.0

    diversity_score = 1 - avg_similarity
    return diversity_score