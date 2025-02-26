import csv
import os

def log_metrics(run_id, 
                timestamp, 
                runtime, 
                hyperparameters, 
                fitness_results, 
                validation_results,
                mem_usage,
                cpu_usage,
                ):
    # connect to metrics log
    cd = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(cd, "metrics_log.csv")
    
    csv_row = {
        # Write out the current run's metadata
        "run_id": run_id,
        "timestamp": timestamp,
        "runtime": runtime,
        
        # Write out hyperparameters from config
        "population_size": hyperparameters["population_size"],
        "num_generations": hyperparameters["num_generations"],
        "mutation_method": hyperparameters["mutation_method"],
        "crossover_method": hyperparameters["crossover_method"],
        #"num_elite_individuals": hyperparameters["num_elite_individuals"],
        #"num_elite_groups": hyperparameters["num_elite_groups"],
        "fitness_function": hyperparameters["fitness_function"],
        #"num_threads": hyperparameters["num_threads"],
        
        # Write out fitness and validation results
        "fitness_results": fitness_results,
        "validation_results": validation_results,
        
        # Write out performance metrics
        "mem_usage": mem_usage,
        "cpu_usage": cpu_usage
    }
    file_exists = os.path.exists(file_path)
    with open(file_path, 'a', newline="") as file:
        writer = csv.DictWriter(file, fieldnames=csv_row.keys())
        if not file_exists:
            writer.writeheader()
        writer.writerow(csv_row)