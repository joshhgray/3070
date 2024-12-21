from src.evolutionary_system.utils.config_loader import load_config
from src.evolutionary_system.utils.logger import log_metrics
from src.evolutionary_system.ga import run_ga
from src.data_pipeline.make_graph import make_graph
from src.data_pipeline.load_population import load_initial_population
import datetime
import psutil # type: ignore
import uuid
import time


# TODO - move into src/evolutionary_systm

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
    initial_population = make_graph(load_initial_population())
    
    final_population = run_ga(
        initial_population=initial_population,
        population_size=config["population_size"],
        num_generations=config["num_generations"],
        mutation_rate=config["mutation_rate"],
        crossover_rate=config["crossover_rate"],
        num_elite_individuals=config["num_elite_individuals"],
        num_elite_groups=config["num_elite_groups"],
        selection_method=config["selection_method"],
        fitness_weights=config["fitness_weights"],
        num_threads=config["num_threads"]
    )
    
    # TODO - append other fitness results
    fitness_results = []
    fitness_results.append(float(final_population.nodes["Population"]["diversity"]))
    
    # End benchmarking
    end_time = time.time()
    runtime = end_time - start_time
    # using oneshot here because I plan on running additional processing benchmarks
    with process.oneshot():
        cpu_usage = process.cpu_percent(interval=None)
        # TODO - implement continuous memory tracking
    
    # TODO - fill rest of parameters
    log_metrics(run_id,
                timestamp,
                runtime=runtime,
                hyperparameters=config,
                fitness_results=fitness_results,
                validation_results={}, # TODO - handle val res
                mem_usage=0, # TODO - handle mem calc
                cpu_usage=cpu_usage,
                )

#def stop_ga():
    # TODO - implement stopping feature


    
#start_ga()