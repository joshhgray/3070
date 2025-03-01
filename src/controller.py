from src.evolutionary_system.mutation_operations.hydroxylate_mutate import hydroxylate_mutate
from src.evolutionary_system.mutation_operations.atomic_substitution import atomic_substitution
from src.evolutionary_system.selection_operations.rank_based_selection import rank_based_selection
from src.evolutionary_system.selection_operations.verhulst_based_selection import verhulst_based_selection
from src.evolutionary_system.selection_operations.stochastic_universal_sampling import stochastic_universal_sampling
from src.evolutionary_system.fitness_operations.calculate_qed import calculate_qed
from src.evolutionary_system.utils.config_loader import load_config
from src.evolutionary_system.utils.logger import log_metrics
from src.evolutionary_system.ga import run_ga
from src.data_pipeline.make_population_graph import make_population_graph
from rdkit import RDLogger
import datetime
import psutil
import uuid
import time

# Suppress RDKit warnings globally
#RDLogger.DisableLog('rdApp.*')

def start_ga():
    """
    Load population, Initialize and Run the Genetic Algorithm (GA) and
    facilitate benchmarking and logging of metrics.
    """

    # Generate primary key run ID
    run_id = str(uuid.uuid4())
    start_time = time.time()
    timestamp = datetime.datetime.now().isoformat()

    # Start tracking CPU usage
    process = psutil.Process()
    process.cpu_percent(interval=None)
    
    # Load GA configuration
    config = load_config()

    # Mutation options mapped
    # TODO - unify naming convention to eliminate map
    mutation_method_map = {
        "hydroxylate": hydroxylate_mutate,
        "atomic_substitution": atomic_substitution
    }

    # TODO Crossover options mapped
    # crossover_method_map = {}

    # Selection options mapped
    selection_method_map = {
        "verhulst": verhulst_based_selection,
        "ranked": rank_based_selection,
        # TODO - Add more
    }

    # Fitness options mapped
    fitness_function_map = {
        "qed": calculate_qed,
    }

    # Assign user selected methods and functions
    selected_mutations = config["mutation_methods"]
    mutation_methods = [mutation_method_map[method] for method in selected_mutations]

    #crossover_method = crossover_method_map[config["crossover_method"]] TODO
    selection_method = selection_method_map[config["selection_method"]]
    fitness_function = fitness_function_map[config["fitness_function"]]

    
    # Load population
    print("Loading Population Graph...")
    population_size = config["population_size"]
    initial_population_graph = make_population_graph(population_size)
    
    
    print("Running GA...")
    final_population, diversity_log = run_ga(
        initial_population=initial_population_graph,
        population_size=population_size,
        num_generations=config["num_generations"],
        mutation_methods=mutation_methods,
        crossover_method=None, # TODO
        num_elite_individuals=None, # TODO - backlog
        num_elite_groups=None, # TODO - backlog
        selection_method=selection_method,
        fitness_function=fitness_function,
        num_threads=None, # TODO - backlog
    )
    
    # Collect fitness results
    # TODO - append other fitness results
    fitness_results = {'diversity': diversity_log}
    
    # End benchmarking
    end_time = time.time()
    runtime = end_time - start_time
    # using oneshot here because I plan on running additional processing benchmarks
    with process.oneshot():
        cpu_usage = process.cpu_percent(interval=None)
        # TODO - implement continuous memory tracking
    
    # Log results
    # TODO - fill rest of parameters
    print("Logging Metrics...")
    log_metrics(run_id,
                timestamp,
                runtime=runtime,
                hyperparameters=config,
                fitness_results=fitness_results,
                validation_results={}, # TODO - handle val res
                mem_usage=0, # TODO - handle mem calc
                cpu_usage=cpu_usage,
                )
    print("Simulation complete.")

    return final_population, diversity_log