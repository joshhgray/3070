import random
import numpy as np
from src.evolutionary_system.selection_operations.selection import Selection

class StochasticUniversalSampling(Selection):


    def select(self, population):
        """
        Performs Stochastic Universal Sampling to select parents.

        :param population:
        :param carrying_capacity:
        :returns
        """
        # Calculate total fitness of population and normalize (prevent divide by zero)
        # Conditional prevents access to group and population level nodes
        individuals = [node for node in population.nodes if population.nodes[node].get("level") == "Individual"]
        
        total_fitness = sum(population.nodes[node].get("raw_fitness") for node in individuals)

        # Prevent division by 0 from being a possibility
        if total_fitness == 0:
            probabilities = np.ones(len(individuals))
        else:
            probabilities = np.array([population.nodes[node]["raw_fitness"] / total_fitness for node in individuals])

        """Roulette-wheel Selection"""
        # Select a minimum of 10 mols as parents
        num_parents = max(10, min(len(individuals) // 5, len(individuals)))
        # Selection points
        step_size = 1 / num_parents
        # Random selection of start point
        start_point = random.uniform(0, step_size)
        cumulative_probabilities = np.cumsum(probabilities)

        selected_parents = []
        i = 0
        for _ in range(num_parents):
            while i < len(cumulative_probabilities) and cumulative_probabilities[i] < start_point:
                i += 1
            selected_parents.append(individuals[i])
            start_point += step_size

        return selected_parents