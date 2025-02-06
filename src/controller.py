from src.evolutionary_system.utils.config_loader import load_config
from src.evolutionary_system.utils.logger import log_metrics
from src.evolutionary_system.ga import run_ga
from src.data_pipeline.make_population_graph import make_population_graph
import datetime
import psutil
import uuid
import time

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
    
    # Load population
    print("Loading Population Graph...")
    population_size = config["population_size"]
    initial_population_graph = make_population_graph(population_size)
    
    
    print("Running GA...")
    final_population, diversity_log = run_ga(
        initial_population=initial_population_graph,
        population_size=population_size,
        num_generations=config["num_generations"],
        mutation_rate=config["mutation_rate"],
        crossover_rate=config["crossover_rate"],
        num_elite_individuals=config["num_elite_individuals"],
        num_elite_groups=config["num_elite_groups"],
        selection_method=config["selection_method"],
        fitness_weights=config["fitness_weights"],
        num_threads=config["num_threads"]
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

#def stop_ga():
    # TODO - implement stopping feature
