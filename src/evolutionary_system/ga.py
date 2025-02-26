from src.evolutionary_system.fitness_operations.calculate_qed import calculate_qed
from src.evolutionary_system.fitness_operations.calculate_population_diversity import calculate_population_diversity
from src.evolutionary_system.selection_operations.rank_based_selection import rank_based_selection
from src.evolutionary_system.mutation_operations.hydroxylate_mutate import hydroxylate_mutate
from src.evolutionary_system.utils.nx_graph_to_mol import nx_graph_to_mol
from src.evolutionary_system.utils.ga_state import update_latest_diversity, update_latest_population, is_ga_active, set_ga_active
from src.data_pipeline.mol_to_graph import mol_to_graph 
import plotly.express as px
from rdkit import Chem
import pandas as pd
import uuid
import csv
import os

latest_diversity_log = []
ga_active = False

def run_ga(initial_population, 
           population_size, 
           num_generations, 
           mutation_method, 
           crossover_method, 
           num_elite_individuals,
           num_elite_groups,
           selection_method,
           fitness_function,
           num_threads,
           selection_cutoff=50):
    
    diversity_log = []
    working_population = initial_population

    cd = os.path.dirname(os.path.abspath(__file__))
    diversity_log_path = os.path.join(cd, 'utils/diversity_log.csv')

    """EVALUATE INITIAL FITNESS"""
    for node in working_population.nodes:
        if working_population.nodes[node]["level"] == "Individual":
            # this will eventually just take an all encompasing fitness function
            # where all metrics are calculated at once. maybe
            compounds = working_population.nodes[node].get("compounds", [])
            mol_graph = compounds[0].get("mol_graph") if compounds else None
            working_population.nodes[node]["raw_fitness"] = fitness_function(mol_graph)
            # will be calculated in fitness.py later - this is temporary

    """
    GENERATIONAL LOOP
    """
    for generation in range(num_generations):
        if not is_ga_active():
            print(f"GA manually stopped.")
            break

        parents = selection_method(working_population, selection_cutoff, population_size)

        """
        MUTATION OPERATION LOGIC
        """
        for parent in parents:
            # Extract mol graph from parent
            # important to add backup logic bc it is easy to
            # create invalid compounds throughout the genetic operations
            compounds = working_population.nodes[parent].get("compounds", [])
            mol_graph = compounds[0].get("mol_graph") if compounds else None
            if mol_graph is None:
                continue

            # Convert mol to rw_mol with nx_graph_to_mol
            rw_mol = nx_graph_to_mol(mol_graph, return_rwmol=True)

            # Apply hydroxylate mutation
            mutated_mol = mutation_method(rw_mol)

            # TODO - add logic to add multiple/different mutation ops

            # Continue on if mutation has failed
            if mutated_mol == rw_mol or mutated_mol is None or mutated_mol.GetNumAtoms() == 0:
                continue
            
            try:
                Chem.SanitizeMol(mutated_mol)
            except Exception as e:
                print(f"Santization failed after mutation: {e}")
                continue

            # add successfully created offspring to offspring list
            # Extract new SMILES string 
            mutated_smiles = Chem.MolToSmiles(mutated_mol)
            
            # Convert back to nx.Graph
            mutated_graph = mol_to_graph(mutated_smiles)

            # Formulate and store offspring to be added to population
            new_id = f"mutated_mol_{uuid.uuid4().hex[:10]}"
            working_population.add_node(new_id, level="Individual", compounds=[])
            working_population.nodes[new_id]["compounds"] = [{"mol_graph": mutated_graph}]
            working_population.nodes[new_id]["qed"] = calculate_qed(mutated_graph)
            working_population.nodes[new_id]["raw_fitness"] = working_population.nodes[new_id]["qed"]

        """
        CROSSOVER OPERATION LOGIC
        """
        
        """
        SELECTION OPERATION LOGIC
        """
        parents = selection_method(working_population, selection_cutoff, population_size)

        """
        CALCULATE DIVERSITY & LOG METRICS
        """
        # TODO - pick a better variable name
        n=1 # I've made this a variable for now, I may add the ability for the user to change it - TODO
        # Calculate population-wide diversity
        if generation % n == 0:
            diversity = calculate_population_diversity(working_population)
            working_population.nodes["Population"]["diversity"] = diversity
            diversity_log.append(diversity)

            with open(diversity_log_path, "a") as f:
                writer = csv.writer(f)
                writer.writerow((generation, diversity))

            
        # Save the latest population and diversity score log in global state
        update_latest_population(working_population)
        update_latest_diversity(diversity_log)
    
    # Update GA state upon completion of Generational Loops
    set_ga_active(False)
    return working_population, diversity_log