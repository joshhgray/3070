import random
import numpy as np
from src.evolutionary_system.selection_operations.selection import Selection

class StochasticUniversalSampling(Selection):


    def select(self, population, root_node="Population", selection_cutoff=0.2):
        """
        Performs Stochastic Universal Sampling to select parents.
        Can perform on set of individuals in full population or set of individuals within 
        a given group. 

        Based off of: https://en.wikipedia.org/wiki/Stochastic_universal_sampling
        Original Source: 
        Baker, J. E. (1987, July). Reducing bias and inefficiency in the selection algorithm. 
        In Proceedings of the second international conference on genetic algorithms (Vol. 206, 
        pp. 14-21).

        :param population: NetworkX graph representign the population.
        :param root_node: A group level node.
        :param selection_cutoff: Percentage of group or population to select.
        :returns: List of selected individual node IDs
        """
        # Group level
        if root_node != "Population":
            individuals = list(population.successors(root_node))
        # Population level
        else:
            individuals = [node for node in population.nodes if population.nodes[node].get("level") == "Individual"]
        
        fitnesses = np.array([population.nodes[node].get("raw_fitness") for node in individuals])
        total_fitness = np.sum(fitnesses)

        # RWS requires normalized probabilities
        if total_fitness == 0: # avoid divide zero
            probabilities = np.ones(len(individuals)) / len(individuals)
        else:
            probabilities = fitnesses / total_fitness

        # SUS setup
        selected_individuals = []
        # minimum of 1 to avoid zero division
        num_selected_individuals = max(1, int(len(individuals) * selection_cutoff))

        # Distance between selection points
        p = 1.0 / num_selected_individuals
        # Random start location
        start = random.uniform(0, p)
        # Selection points
        pointers = [start + i * p for i in range(num_selected_individuals)]
        cumulative_probabilities = np.cumsum(probabilities)
        i = 0

        """
        Roulette-wheel Selection
        """
        for pointer in pointers:
            while i < len(cumulative_probabilities) and cumulative_probabilities[i] < pointer:
                i += 1
            if i < len(individuals): # avoid out of bounds
                selected_individuals.append(individuals[i])

        return selected_individuals