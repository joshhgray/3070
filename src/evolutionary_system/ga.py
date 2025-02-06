import random
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

    #print("evaluating initial fitness")
    # Evaluate  Initial Fitness
    for node in working_population.nodes:
        if working_population.nodes[node]["level"] == "Individual":
            # this will eventually just take an all encompasing fitness function
            # where all metrics are calculated at once. maybe
            compounds = working_population.nodes[node].get("compounds", [])
            structure = compounds[0].get("mol_graph") if compounds else None
            working_population.nodes[node]["qed"] = calculate_qed(structure)
            # will be calculated in fitness.py later - this is temporary
            working_population.nodes[node]['raw_fitness'] = working_population.nodes[node]["qed"]

    #print("evaluating fitness")
    for generation in range(num_generations):

        # Selection TODO - apply a proper selection technique
        # Currently just selected top half of population based on raw fitness
        mols = [node for node in working_population.nodes if working_population.nodes[node]["level"] == "Individual"]
        sorted_mols = sorted(mols,
                                    key=lambda x: working_population.nodes[x]["raw_fitness"], 
                                    reverse=True)
        parents = sorted_mols[: len(sorted_mols) // 2]

        # Crossover 
        # TODO make a working crossover function
        new_population = [] 
        for i in range(0, len(parents) -1, 2): # 2 parents
            p1 = working_population.nodes[parents[i]].get("compounds", [])
            p2 = working_population.nodes[parents[i+1]].get("compounds", [])
            # TODO - temporarily this takes both p1 an p1 effectively making it do nothing (since it's not implemented properly yet)
            child1, child2 = p1, p2 #crossover(p1, p1)
            new_population.append(child1)
            new_population.append(child2)

        # Update graph (Build new population)
        for i, mol in enumerate(new_population):
            new_id = f"evoMol_{uuid.uuid4().hex[:10]}"
            working_population.add_node(new_id, level="Individual", compounds=[])
            working_population.nodes[new_id]["comopounds"] = mol
            working_population.nodes[new_id]["qed"] = calculate_qed(mol)
            working_population.nodes[new_id]['raw_fitness'] = working_population.nodes[new_id]["qed"]
        
        
        # Calculate population-wide diversity
        diversity = calculate_diversity(working_population)
        working_population.nodes["Population"]["diversity"] = diversity
        diversity_log.append(diversity)
    
    return working_population, diversity_log