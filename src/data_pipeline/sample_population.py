import random

def sample_bgcs(bgc_data, population_size):
    """
    Randomly samples a population of the given size from the stored 
    pool of BGCs.

    :param bgc_data: List of dictionaries - each holding the data and metadata of
                     a BGC in the population.
    :param population_size: The desired size of the population to be sent to the GA.
    """
    return random.sample(bgc_data, population_size)