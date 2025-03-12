from src.evolutionary_system.fitness_operations.aggregate_fitness import aggregate_fitness
from src.evolutionary_system.fitness_operations.calculate_population_diversity import calculate_population_diversity
from src.evolutionary_system.selection_operations.rank_based_selection import rank_based_selection
from src.evolutionary_system.mutation_operations.hydroxylate_mutate import hydroxylate_mutate
from src.evolutionary_system.utils.nx_graph_to_mol import nx_graph_to_mol
from src.evolutionary_system.utils.ga_state import (
    update_latest_diversity, update_latest_population, is_ga_active, 
    set_ga_active, update_latest_crossover_rates, update_latest_mutation_rates, 
    update_crossover_log, update_mutation_log, update_current_population_size, 
    update_current_generation_number
)
from src.data_pipeline.mol_to_graph import mol_to_graph
import plotly.express as px
from rdkit import Chem
import pandas as pd
import random
import uuid
import csv
import os

# Initialize Logs
latest_diversity_log = []

ga_active = False

# Probabilities of any of the mutations to occur in current generational cycle
# TODO - eventually get this to front-end 
MUTATION_PROBABILITIES = {
    "hydroxylate_mutate": 0.2,
    "methylate_mutate": 0.1,
    "atomic_substitution": 0.4,
}

def apply_mutation(working_population, mutation_methods, parents):
    """
    Applies probabalistic mutation to the parent molecules.
    Continusously tracks success rate of mutation operations.

    :param working_population: Current population graph (nx.DiGraph)
    :param mutation_methods: List of possible mutation methods 
    :param parents: List of selected parent nodes
    :returns: List of mutation success and attempt rates.
    """
    mutation_attempts = 0
    mutation_successes = 0
    for parent in parents:
        # Extract mol graph from parent
        # important to add backup logic bc it is easy to
        # create invalid compounds throughout the genetic operations
        compounds = working_population.nodes[parent].get("compounds")
        mol_graph = compounds[0].get("mol_graph") if compounds else None
        if mol_graph is None:
            continue
        
        # Convert mol to rw_mol with nx_graph_to_mol
        rw_mol = nx_graph_to_mol(mol_graph, return_rwmol=True)

        if rw_mol is None:
            continue

        # mutated_mol must be initialized pre-loop in the event a molecule gets no mutation
        mutated_mol = rw_mol

        # Apply each mutation method in sequence
        for mutation_method in mutation_methods:
            # Retrieve mutation method by name (string)
            if random.random() < MUTATION_PROBABILITIES[mutation_method.__name__]:
                mutation_attempts += 1
                mutated_mol = mutation_method(rw_mol)

                # Revert to original if mutation or sanitization fails
                if mutated_mol is None or mutated_mol.GetNumAtoms() == 0:
                    mutated_mol = rw_mol
                try:
                    Chem.SanitizeMol(mutated_mol)
                except Exception as e:
                    mutated_mol = rw_mol

                # Check if mutation was successful
                mutated_smiles = Chem.MolToSmiles(mutated_mol)
                original_smiles = Chem.MolToSmiles(rw_mol)
                if mutated_smiles != original_smiles:
                    mutation_successes += 1
        
        # New mol data
        mutated_smiles = Chem.MolToSmiles(mutated_mol)
        mutated_graph = mol_to_graph(mutated_smiles)
        if mutated_graph is None:
            continue

        # TODO - would it be better to make an add_to_population function that handles this stuff?
        # Formulate and store offspring to be added to population
        new_id = f"offspring_{uuid.uuid4().hex[:10]}"
        working_population.add_node(new_id, level="Individual", compounds=[])
        working_population.nodes[new_id]["compounds"] = [{"structure": mutated_smiles, "mol_graph": mutated_graph}]
        working_population.nodes[new_id]["raw_fitness"] = aggregate_fitness(mutated_graph)
    
    return mutation_attempts, mutation_successes

def apply_crossover(working_population, crossover_methods, parents):
    """
    Pairs parents together and iteratively attempts to perform genetic crossover, creating new
    offspring. If validation of molecule is successful, the offspring is added to the population.
    Continuosly tracks success rate of crossover mutations.

    :param working_population: Current Population graph (nx.DiGraph)
    :param crossover_methods: list of possible crossover methods
    :param parents: list of selected parent nodes
    :returns: List of crossover attempt and success rates.
    """

    # Parent pairing logic 
    # TODO - Once grouping is implemented update logic, for now pairing logic is just random
    random.shuffle(parents)
    # Create pairs of parents
    parent_pairs = [(parents[i], parents[i+1]) for i in range(0, len(parents) - 1, 2)]

    crossover_attempts = 0
    crossover_successes = 0

    # Iterate over parents
    for parent1, parent2 in parent_pairs:
        # Extract each parents molecular graph
        mol_graph1 = working_population.nodes[parent1].get("compounds")[0].get("mol_graph")
        mol_graph2 = working_population.nodes[parent2].get("compounds")[0].get("mol_graph")

        crossover_attempts += 1
        # Apply Crossover technique
        # TODO - currently just doing a random choice, want to make this probabalistic like
        #        mutation, gotta figure a few things out first tho. (also there is only 1 x-over op rn)
        crossover_method = random.choice(crossover_methods)
        offspring_graph = crossover_method(mol_graph1, mol_graph2)

        # Skip if offspring is invalid
        if not offspring_graph:
            continue

        # Convert to regular RDKit mol
        offspring_mol = nx_graph_to_mol(offspring_graph, return_rwmol=False)

        # Confirm chemical validity, Skip if offspring is invalid
        try:
            Chem.SanitizeMol(offspring_mol)
        except Exception:
            continue

        
        crossover_successes += 1
        # Store new molecule in population
        offspring_smiles = Chem.MolToSmiles(offspring_mol)
        offspring_graph_valid = mol_to_graph(offspring_smiles)
        if offspring_graph_valid:
            new_id = f"offspring_{uuid.uuid4().hex[:10]}"
            working_population.add_node(new_id, level="Individual", compounds=[])
            working_population.nodes[new_id]["compounds"] = [{"structure": offspring_smiles, "mol_graph": offspring_graph_valid}]
            # TODO - move this to GA so that fitness can be properly evaluated for new mol (or do here idk)
            working_population.nodes[new_id]["raw_fitness"] = aggregate_fitness(offspring_graph_valid)

    return crossover_attempts, crossover_successes

def apply_selection(working_population, selection_method, selection_cutoff, carrying_capacity):
    """
    Applies the selected selection method to choose parents for reproduction, assigning corresponding
    parameters.
    
    :param working_population:
    :param selection_method:
    :param selection_cutoff:
    :param carrying_capacity:
    :returns: List of selected molecules
    """
    if selection_method.__name__ == "stochastic_universal_sampling":
        parents = selection_method(working_population, selection_cutoff, carrying_capacity)
    else:
        parents = selection_method(working_population, selection_cutoff)

    return parents

def run_ga(initial_population, population_size, num_generations, 
           mutation_methods, crossover_methods, num_elite_individuals,
           num_elite_groups, selection_method, num_threads, 
           selection_cutoff, carrying_capacity):
    
    # Initialize
    diversity_log = []
    current_generation_number = 0
    current_population_size = population_size
    working_population = initial_population

    cd = os.path.dirname(os.path.abspath(__file__))
    diversity_log_path = os.path.join(cd, 'utils/diversity_log.csv')

    """EVALUATE INITIAL FITNESS"""
    for node in working_population.nodes:
        if working_population.nodes[node]["level"] == "Individual":
            compounds = working_population.nodes[node].get("compounds")
            mol_graph = compounds[0].get("mol_graph") if compounds else None
            working_population.nodes[node]["raw_fitness"] = aggregate_fitness(mol_graph)

    """
    GENERATIONAL LOOP
    """
    for generation in range(num_generations):
        if not is_ga_active():
            break
        current_generation_number += 1

        """
        SELECTION OPERATION 
        """
        parents = apply_selection(working_population, selection_method, selection_cutoff, carrying_capacity)

        """
        MUTATION OPERATION
        """
        mutation_attempts, mutation_successess = apply_mutation(working_population, mutation_methods, parents)

        """
        CROSSOVER OPERATION 
        """
        crossover_attempts, crossover_successess = apply_crossover(working_population, crossover_methods, parents)

        # Log generational mutation and crossover success rates
        if mutation_attempts > 0:
            mutation_rate = (mutation_successess / mutation_attempts * 100)
        if crossover_attempts > 0:
            crossover_rate = (crossover_successess / crossover_attempts * 100)

        # Log to global state TODO - can probably condense these
        update_mutation_log(mutation_rate)
        update_crossover_log(crossover_rate)

        update_latest_mutation_rates(mutation_rate)
        update_latest_crossover_rates(crossover_rate)
        
        """
        CALCULATE DIVERSITY & LOG METRICS
        """
        diversity = calculate_population_diversity(working_population)
        working_population.nodes["Population"]["diversity"] = diversity
        diversity_log.append(diversity)

        with open(diversity_log_path, "a") as f:
            writer = csv.writer(f)
            writer.writerow((generation, diversity))

            
        # Save the latest population and diversity score log in global state
        update_latest_population(working_population)
        update_latest_diversity(diversity_log)
        current_population_size = sum(1 for node in working_population.nodes if working_population.nodes[node]["level"] == "Individual")

        update_current_generation_number(current_generation_number)
        update_current_population_size(current_population_size)
    
    # Update GA state upon completion of Generational Loops
    set_ga_active(False)
    return working_population, diversity_log