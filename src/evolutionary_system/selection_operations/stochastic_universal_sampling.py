import random

def stochastic_universal_sampling(population, population_size):
    """
    Performs Stochastic Universal Sampling to select parents
    """
    bgcs = [node for node in population.nodes if population.nodes[node]["level"] == "Individual"]
    fitness_scores = [population.nodes[node].get("raw_fitness", 0.0) for node in bgcs]

    total_fitness = sum(fitness_scores)
    #TODO