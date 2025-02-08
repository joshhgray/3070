from src.evolutionary_system.fitness_operations.calculate_qed import calculate_qed
from src.evolutionary_system.fitness_operations.calculate_population_diversity import calculate_population_diversity
from src.evolutionary_system.selection_operations.rank_based_selection import rank_based_selection
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
            structure = compounds[0].get("mol_graph") if compounds else None
            working_population.nodes[node]["qed"] = calculate_qed(structure)
            # will be calculated in fitness.py later - this is temporary
            working_population.nodes[node]['raw_fitness'] = working_population.nodes[node]["qed"]

    """
    GENERATIONAL LOOP
    """
    for generation in range(num_generations):

        """
        SELECTION
        """
        parents = rank_based_selection(working_population, selection_cutoff)

        """
        CROSSOVER
        """

        # TODO make a working crossover function
        new_population = [] 
        for i in range(0, len(parents) -1, 2): # 2 parents
            # Extract parent's molecular graphs
            p1 = working_population.nodes[parents[i]].get("compounds", [])
            p2 = working_population.nodes[parents[i+1]].get("compounds", [])
            p1_mol_graph = p1[0].get("mol_graph") if p1 else None
            p2_mol_graph = p2[0].get("mol_graph") if p2 else None

            # TODO - This currently does nothing - just here to test integration
            child1 = p1_mol_graph
            child2 = p2_mol_graph

            new_population.append(child1)
            new_population.append(child2)

        """
        UPDATE GRAPH - BUILD NEW POPULATION
        """
        for i, mol in enumerate(new_population):
            new_id = f"evoMol_{uuid.uuid4().hex[:10]}"
            working_population.add_node(new_id, level="Individual", compounds=[])

            if mol is None:
                print(f"no mol_graph for {new_id}")

            working_population.nodes[new_id]["compounds"] = [{"mol_graph": mol}]
            working_population.nodes[new_id]["qed"] = calculate_qed(mol)
            working_population.nodes[new_id]['raw_fitness'] = working_population.nodes[new_id]["qed"]
        
        """
        CALCULATE DIVERSITY
        """
        # Calculate population-wide diversity
        diversity = calculate_population_diversity(working_population)
        working_population.nodes["Population"]["diversity"] = diversity
        diversity_log.append(diversity)
    
    return working_population, diversity_log