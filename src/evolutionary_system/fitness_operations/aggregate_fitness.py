from src.evolutionary_system.fitness_operations.fitness_functions import (
    calculate_qed, calculate_sa_score, calculate_molecular_weight,
    lipinski_score
)
from src.evolutionary_system.utils.ga_state import get_mol_weight_threshold, get_selected_fitness_functions
import numpy as np
from scipy.stats import hmean

def aggregate_fitness(compound):
    """
    Calculates an aggregate fitness scores based on user selected weights
    for QED, SA Score, Ro5, and Molecular Weight.

    :param compound: Netowrkx Graph Object representing molecular compound.
    :param weights: Dictionary containing user-selected weights for each fitness function.
    :returns: Aggregated fitness score (float between 0.0-1.0).
    """
    # Default to equally weighted aggregate function

    selected_fitness_functions = get_selected_fitness_functions()

    fitness_scores = []

    # QED
    if "qed" in selected_fitness_functions:
        fitness_scores.append(calculate_qed(compound))

    # SA - Score inversed since higher is worse
    if "sa" in selected_fitness_functions:
        sa_score = calculate_sa_score(compound)
        sa_score_inversed_normalized = 1 - (sa_score / 10)
        fitness_scores.append(sa_score_inversed_normalized)
    
    if "ro5" in selected_fitness_functions:
        fitness_scores.append(lipinski_score(compound))


    # Smooth molecular weight threshold penalty w/ Sigmoid Function
    mol_weight = calculate_molecular_weight(compound)
    threshold = get_mol_weight_threshold()
    mol_weight_penalty = 1 / (1 + np.exp(1 + np.exp(0.01 * (mol_weight - threshold))))

    if not fitness_scores:
        return 0.0

    # Maintain the option to select just a single fitness function
    if len(fitness_scores) == 1:
        return round(fitness_scores[0] * mol_weight_penalty, 4)

    mean_fitness_score = hmean(fitness_scores)

    return round(mean_fitness_score * mol_weight_penalty, 4)

