import random
import numpy as np

def verhulst_population_control(population_size, carrying_capacity):
    """
    Control population size using Verhulst's model of population growth
    dN/dt = rN(1-N/K)
    
    :param population_size: Total number of molecules in population (N)
    :param carrying_capacity: Max population size (K) - 
    TODO - need to find a K that isn't too computationally expensive but still allows exploratory growth 
    :returns: float from 0-1 representing survival probability for mols in current generation
              (aka scale of population growth for current generation)
    """
    # TODO - find good scales - trying out an explicit limit
    if population_size >= carrying_capacity:
        return 0
    else:
        return max(0, 1 - (population_size / (2 * carrying_capacity)))

def stochastic_universal_sampling(population, num_selected, carrying_capacity):
    """
    Performs Stochastic Universal Sampling to select parents 
    with Verhulst/Logistic-based population control

    :param population:
    :param num_selected: Number of molecules to select
    :param carrying_capacity:
    :returns
    """
    # Calculate total fitness of population and normalize (prevent divide by zero)
    # Conditional prevents access to group and population level nodes
    mols = [node for node in population.nodes if population.nodes[node].get("level") == "Individual"]

    total_fitness = sum(population.nodes[node].get("raw_fitness") for node in mols)

    # Prevent division by 0 from being a possibility
    if total_fitness == 0:
        selected_mols = random.sample(mols, min(num_selected), len(mols))
        unselected_mols = set(mols) - set(selected_mols)
        if unselected_mols:
            population.remove_nodes_from(unselected_mols)
        return selected_mols

    # Normalized fitness for each mol
    norm_fits = [population.nodes[node]["raw_fitness"] / total_fitness for node in mols]

    # Use Verhulst population growth function to determine probablity of survival
    survival_propbability = verhulst_population_control(len(population), carrying_capacity)
    if survival_propbability == 0:
        return []

    # Probabilities
    adjusted_probs = [fitness * survival_propbability for fitness in norm_fits]

    total_probs = sum(adjusted_probs)
    if total_probs == 0:
        return []
    
    norm_probs = [prob / total_probs for prob in adjusted_probs]
    cumulative_probabilities = np.cumsum(norm_probs)

    """Roulette-wheel Selection"""
    # Selection points
    step_size = 1 / num_selected 

    # Random selection of start point
    start_point = random.uniform(0, step_size)

    # Pointers
    pointers = [start_point + i * step_size for i in range(num_selected)]

    # *Spin the Wheel*
    selected_mols = []
    i = 0
    for pointer in pointers:
        while i < len(cumulative_probabilities) - 1 and cumulative_probabilities[i] < pointer:
            i += 1
        selected_mols.append(mols[i])
    
    # Some invalids were making it through somehow - TODO - find out why 
    valid_parents = [node for node in selected_mols if population.nodes[node].get("compounds")]

    # Pruning of unselected Inidividuals (Required to prevent exponential growth of population)
    unselected_mols = set(mols) - set(valid_parents)
    if unselected_mols:
        population.remove_nodes_from(unselected_mols)

    return valid_parents
        
    



    