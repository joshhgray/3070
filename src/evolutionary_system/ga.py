import random
from src.evolutionary_system.fitness import calculate_qed, calculate_diversity

def mutate(smiles_str):
    # TODO: Basic mutation
    return smiles_str

def crossover(p1, p2):
    # Basic crossover operation, splitting at random location and combining halves
    crossover_point = random.randint(1, len(p1) - 1)
    child1 = p1[:crossover_point] + p2[crossover_point:]
    child2 = p2[:crossover_point] + p1[crossover_point:]
    return child1, child2 

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
    global is_running
    diversity_log = []

    for generation in range(num_generations):

        working_population = initial_population
        # Evaluate Fitness
        for node in working_population.nodes:
            # Only need to access individuals here
            if working_population.nodes[node]["level"] == "Individual":
                smiles_str = working_population.nodes[node]['smiles_str']
                # this will eventually just take an all encompasing fitness function
                # where all metrics are calculated at once
                working_population.nodes[node]["qed"] = calculate_qed(smiles_str)
                # will be calculated in fitness.py later - this is temporary
                working_population.nodes[node]['raw_fitness'] = working_population.nodes[node]["qed"]

        # top half of population
        individuals = [node for node in working_population.nodes if working_population.nodes[node]["level"] == "Individual"]
        sorted_individuals = sorted(individuals,
                                    key=lambda x: working_population.nodes[x]["raw_fitness"], 
                                    reverse=True)
        parents = sorted_individuals[: len(sorted_individuals) // 2]

        #crossover operation
        new_population = [] 
        for i in range(0, len(parents) -1, 2): # 2 parents
            p1 = working_population.nodes[parents[i]]["smiles_str"]
            p2 = working_population.nodes[parents[i+1]]["smiles_str"]
            child1, child2 = crossover(p1, p1)
            new_population.append(child1)
            new_population.append(child2)

        # graph update, build new pop
        for i, individual in enumerate(new_population):
            node_name = f"Individual_1_{i + 1}"
            if node_name in working_population.nodes:
                working_population.nodes[node_name]["smiles_str"] = individual
        
        working_population.nodes[node]
        
        # Calculate population-wide diversity
        diversity = calculate_diversity(working_population)
        working_population.nodes["Population"]["diversity"] = diversity

        diversity_log.append(diversity)
        
    return working_population, diversity_log

