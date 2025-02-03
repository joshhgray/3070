import random
import pickle
import os

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