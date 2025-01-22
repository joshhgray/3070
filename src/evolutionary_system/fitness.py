# Base fitness function
from rdkit import Chem #type: ignore
from rdkit.Chem import QED, AllChem #type: ignore
from rdkit.DataStructs import TanimotoSimilarity #type: ignore
import numpy as np #type: ignore

def calculate_qed(structure):
    """
    Calculate QED score of a given SMILES string

    :param structure: the SMILES-notation chemical notation of the molecule
    :returns: QED score of the given structure
    """
    try:
        mol = Chem.MolFromSmiles(structure)
        if mol: # Valid compound
            return QED.qed(mol)
        else:
            print(f"Invalid compound: {structure}.")
            return 0.0

    except Exception as e:
        print(f"Error calculating QED for {structure}: {e}")
        return 0.0











# TODO - Archive 
# def calculate_diversity(population):
#     # TODO : Satisfy deprecation warning --> [20:43:06] DEPRECATION WARNING: please use MorganGenerator
#     # Extract list of molecules from smiles strings
#     individuals = [
#         Chem.MolFromSmiles(population.nodes[node]["smiles_str"])
#         for node in population.nodes
#         if "smiles_str" in population.nodes[node]
#         and population.nodes[node]["smiles_str"] is not None
#     ]
    
#     # Generage fingerprints
#     fingerprints = [AllChem.GetMorganFingerprintAsBitVect(individual, 2, nBits=1024) for individual in individuals]
    
#     # Calculate Tanimoto Similarity (Jacdard Index)
#     similarity_scores = []
#     for i in range(len(fingerprints)):
#         for j in range(i+1, len(fingerprints)):
#             similarity_score = TanimotoSimilarity(fingerprints[i], fingerprints[j])
#             similarity_scores.append(similarity_score)
            
#     avg_similarity = float(np.mean(similarity_scores)) if similarity_scores else 0
#     population_diversity = 1 - avg_similarity
#     return population_diversity