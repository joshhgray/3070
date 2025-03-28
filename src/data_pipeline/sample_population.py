import random
import pickle
import os

# TODO - could probably combine these into one class
# TODO - expand to include SELFIES strings
# TODO - expand to include other commonly used forms of chemical representations

def sample_bgcs(population_size):
    """
    Randomly samples a population of the given size from the stored 
    pool of BGCs.

    :param bgc_data: file path to preprocessed BGC data
    :param population_size: The desired size of the population to be sent to the GA.
    :returns sampled_bgcs: 
    """
    
    cd = os.path.dirname(os.path.abspath(__file__))
    population_pool = os.path.join(cd, "../../preprocessed_bgcs.pkl")
    with open(population_pool, "rb") as f:
        bgc_data = pickle.load(f)

    if population_size > len(bgc_data):
        population_size = len(bgc_data)

    sampled_bgcs = random.sample(bgc_data, population_size)

    return sampled_bgcs

def sample_smiles(smiles_list, population_size):
    """
    Randomly samples a population of the given size form a list of SMILES strings.
    
    :param smiles_list: Python list of SMILES strings.
    :param population_size: User-selected population size.
    :returns: List of randomly sampled SMILES matching the size of the population.
    """
    if population_size > len(smiles_list):
        population_size = len(smiles_list)
    
    sampled_smiles = random.sample(smiles_list, population_size)

    return sampled_smiles