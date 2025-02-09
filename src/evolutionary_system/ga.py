from src.evolutionary_system.fitness_operations.calculate_qed import calculate_qed
from src.evolutionary_system.fitness_operations.calculate_population_diversity import calculate_population_diversity
from src.evolutionary_system.selection_operations.rank_based_selection import rank_based_selection
from src.evolutionary_system.mutation_operations.hydroxylate_mutate import hydroxylate_mutate
from src.evolutionary_system.utils.nx_graph_to_mol import nx_graph_to_mol
from src.data_pipeline.mol_to_graph import mol_to_graph 
from rdkit import Chem
import uuid

def run_ga(initial_population, 
           population_size, 
           num_generations, 
           mutation_rate, 
           crossover_rate, 
           num_elite_individuals,
           num_elite_groups,
           selection_method,
           fitness_weights,
           num_threads,
           selection_cutoff=50):
    
    diversity_log = []
    working_population = initial_population

    """EVALUATE INITIAL FITNESS"""
    for node in working_population.nodes:
        if working_population.nodes[node]["level"] == "Individual":
            # this will eventually just take an all encompasing fitness function
            # where all metrics are calculated at once. maybe
            compounds = working_population.nodes[node].get("compounds", [])
            mol_graph = compounds[0].get("mol_graph") if compounds else None
            working_population.nodes[node]["qed"] = calculate_qed(mol_graph)
            # will be calculated in fitness.py later - this is temporary
            working_population.nodes[node]['raw_fitness'] = working_population.nodes[node]["qed"]

    """
    GENERATIONAL LOOP
    """
    for generation in range(num_generations):

        """
        SELECTION OPERATION LOGIC
        """
        parents = rank_based_selection(working_population, selection_cutoff)

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
            mutated_mol = hydroxylate_mutate(rw_mol)

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
        TRUNCATE POPULATION 
        """
        while len(working_population.nodes) > population_size:
            # Delete worst-performing individuals based, just based on raw fitness for now - TODO
            worst_mol = min(working_population.nodes, key=lambda x: working_population.nodes[x].get('raw_fitness', float("inf")))
            working_population.remove_node(worst_mol)

        """
        CALCULATE DIVERSITY
        """
        # Calculate population-wide diversity
        diversity = calculate_population_diversity(working_population)
        working_population.nodes["Population"]["diversity"] = diversity
        diversity_log.append(diversity)
    
    return working_population, diversity_log