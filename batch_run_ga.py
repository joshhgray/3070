from tqdm.contrib.concurrent import process_map
from multiprocessing import Pool, cpu_count
from rdkit import RDLogger
import itertools
import uuid

# Supress RDKit logging
RDLogger.DisableLog('rdApp.*')

"""
To run: python3 batch_run_ga.py | tee batch_run.log
"""

"""
DEFINE PARAMETERS
"""
# carrying_capacities = [2000, 5000]
# crossover_methods = [["mcs_crossover"]]
# data_sources = ["smiles", "bgc"]
# elitism_weights = [1, 2]
# mutation_weights = [
#     {"atomic_substitution": 0.25, "functional_group": 0.25, "ring": 0.25, "fragment": 0.25},
#     {"atomic_substitution": 0.75, "functional_group": 0.25, "ring": 0.25, "fragment": 0.25},
#     {"atomic_substitution": 0.25, "functional_group": 0.75, "ring": 0.25, "fragment": 0.25},
#     {"atomic_substitution": 0.25, "functional_group": 0.25, "ring": 0.75, "fragment": 0.25},
#     {"atomic_substitution": 0.25, "functional_group": 0.25, "ring": 0.25, "fragment": 0.75},
# ]
# num_generations = [1000, 2000]
# population_sizes = [1000, 2000, 5000]
# selection_methods = ["stochastic_universal_sampling", "rank_based_selection"]

"""
DEMO PARAMETERS - quick run
"""
carrying_capacities = [2000]
crossover_methods = [["mcs_crossover"]]
data_sources = ["smiles"]
elitism_weights = [1]
mutation_weights = [
    {"atomic_substitution": 0.25, "functional_group": 0.25, "ring": 0.25, "fragment": 0.25},
]
num_generations = [100, 1000]
population_sizes = [1000]
selection_methods = ["stochastic_universal_sampling"]

"""
GENERATE CONFIGS
"""
def generate_configs():
    configs = []
    for (capacity, crossover, source, elitism, mutation,
         # compile product of each configuration
         generation, pop_size, selection) in itertools.product(
             carrying_capacities, crossover_methods, data_sources,
             elitism_weights, mutation_weights, num_generations,
             population_sizes, selection_methods):
        
        # Filter out invalid configurations
        if pop_size > capacity:
            continue
        
        # build config
        config = {
            "run_id": str(uuid.uuid4())[:10], # add custom id for each run
            "carrying_capacity": capacity, 
            "crossover_methods": crossover,
            "data_source": source,
            "elitism_weight": elitism,
            "mutation_weights": mutation,
            "num_generations": generation,
            "population_size": pop_size,
            "selection_method": selection,
            }
        configs.append(config)
    return configs

""" 
RUN SIMULATIONS
"""
def run_ga(config):
    #print(f"Running GA {config['run_id']} with a population size of {config['population_size']} and {config['num_generations']} generations.")
    try:
        from src.evolutionary_system.GeneticAlgorithm import GeneticAlgorithm
        from src.evolutionary_system.utils.ga_state import set_ga_active
        set_ga_active(True)
        ga = GeneticAlgorithm(config)
        ga.start_ga()
    except Exception as e:
        print(f"ERROR in run: {e}")

import multiprocessing as mp
if __name__ == "__main__": 
    mp.set_start_method("spawn", force=True)
    configs = generate_configs()
    num_threads = min(cpu_count(), len(configs))# run on max available threads for maximum efficiency
    
    print(f"Running {len(configs)} batch jobs on {num_threads} threads...")

    # tqdm process_map facilitates tracking over multiiprocessing
    process_map(run_ga, configs, max_workers=9, chunksize=1)
    # with Pool(processes=num_threads) as pool:
    #     results = pool.map(run_ga, configs)

    print("Completed batch runs.")
    