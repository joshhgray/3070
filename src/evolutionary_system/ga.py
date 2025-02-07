from src.evolutionary_system.fitness import calculate_qed, calculate_diversity
from rdkit import Chem #type:ignore
import uuid

def mutate():
    pass

def crossover():
    pass

def run_ga(initial_population, 
           population_size, 
           num_generations, 
           mutation_rate, 
           crossover_rate, 
           num_elite_individuals,
           num_elite_groups,
           selection_method,
           fitness_weights,
           num_threads):
    
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
        # TODO - apply a proper selection technique
        # Currently just selected top half of population based on raw fitness
        mols = [node for node in working_population.nodes if working_population.nodes[node]["level"] == "Individual"]
        sorted_mols = sorted(mols,
                                    key=lambda x: working_population.nodes[x]["raw_fitness"], 
                                    reverse=True)
        parents = sorted_mols[: len(sorted_mols) // 2]

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
        diversity = calculate_diversity(working_population)
        working_population.nodes["Population"]["diversity"] = diversity
        diversity_log.append(diversity)
    
    return working_population, diversity_log