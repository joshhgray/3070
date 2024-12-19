from src.evolutionary_system.utils.config_loader import load_config
from src.evolutionary_system.utils.logger import log_metrics
from src.evolutionary_system.ga import run_ga
from src.data_pipeline.make_graph import make_graph
import datetime
import psutil
import uuid
import time

# not sure if i'm going to keep this in here
# will do the same for validation results
def unpack_fitness_results(final_population):
    fitness_results = []
    for node in final_population.nodes:
        if final_population.nodes[node]["level"] == "Individual":
            fitness_results.append({
                "smiles_str": final_population.nodes[node]["smiles_str"],
                "raw_fitness": final_population.nodes[node]["raw_fitness"]
                # TODO - add remaining results
            })    
    
def start_ga():
    # Setup benchmarking tools
    
    # Generate run_id (primary key) with uuid4
    run_id = str(uuid.uuid4())
    start_time = time.time()
    timestamp = datetime.datetime.now().isoformat()
    '''
    Start tracking CPU usage -
    Measured as % of total available CPU from here
    to end of GA execution
    '''
    process = psutil.Process()
    process.cpu_percent(interval=None)
    
    # Set up and run GA
    config = load_config()
    
    population = make_graph()
    final_population = run_ga(
        population=population,
        population_size=config["population_size"],
        num_generations=config["num_generations"],
        mutation_rate=config["mutation_rate"],
        crossover_rate=config["crossover_rate"],
        num_elite_individuals=config["num_elite_individuals"],
        num_elite_groups=config["num_elite_groups"],
        fitness_weights=config["fitness_weights"],
        num_threads=config["num_threads"]
    )
    fitness_results = unpack_fitness_results(final_population)
    
    # End benchmarking
    end_time = time.time()
    runtime = end_time - start_time
    with process.oneshot():
        cpu_usage = process.cpu_percent(interval=None)
        # TODO - implement continuous memory tracking
    
    # TODO - fill rest of parameters
    log_metrics(run_id,
                timestamp,
                runtime=runtime,
                hyperparameters=config,
                fintess_resutls=fitness_results,
                validation_results={}, # TODO - handle val res
                mem_usage=0, # TODO - handle mem calc
                cpu_usage=cpu_usage,
                )