from src.evolutionary_system.fitness_operations.fitness_functions import (
    calculate_qed, calculate_sa_score, calculate_molecular_weight,
    lipinski_score, calculate_logp,
)
from src.evolutionary_system.utils.ga_state import (
    get_mw_range, get_mw_target, get_active_filters, 
    get_tuning_weights, get_logp_range
)
import numpy as np
from scipy.stats import hmean

def aggregate_fitness(mol):
    """
    Calculates an aggregate fitness scores based on user selected weights
    for QED, SA Score, Ro5, and Molecular Weight.

    :param mol: RDKit Mol Object.
    :param weights: Dictionary containing user-selected weights for each fitness function.
    :returns: Aggregated fitness score (float between 0.0-1.0).
    """
    # Default to equally weighted aggregate function
    if mol is None:
        return 0.0

    # Retreive live fitness filters and weights
    filters = get_active_filters()
    weights = get_tuning_weights()

    """
    Filters - hard filters that ensure failing molecules are removed from the population.
    """ 
    if filters is not []:
        # Lipinski's Rule of Five
        if "ro5" in filters and lipinski_score(mol) >= 0.6:
            print("filtered by ro5")
            return 0.01

        # Molecular Weight Filter
        mol_weight = calculate_molecular_weight(mol)
        mw_range = get_mw_range()
        if mol_weight < mw_range[0] or mol_weight > mw_range[1]:
            return 0.01

        # logP Filter
        logp = calculate_logp(mol)
        logp_range = get_logp_range()
        if logp < logp_range[0] or logp > logp_range[1]:
            return 0.01

    """
    Tuning of Weights
    """
    tuning_functions = {
        "qed": lambda mol: calculate_qed(mol),
        "sa": lambda mol: 1 - (calculate_sa_score(mol) / 10), # Inverted and Normalized
        "mol_weight": lambda mol: 1 - min(abs(calculate_molecular_weight(mol) - get_mw_target()) / 250, 1)
    }

    # Combine
    scores = []
    for key, weight in weights.items():
        if weight > 0 and key in tuning_functions:
            score = tuning_functions[key](mol)
            scores.append(score * weight)

    if not scores:
        return 0.0
    
    # Return normalized score
    return sum(scores) / sum(weights.values())



