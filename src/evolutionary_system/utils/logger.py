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
                population_type,
                fitness_config
                ):
    # connect to metrics log
    cd = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(cd, "metrics_log.csv")

    # Extract logs - convert to string for csv compatability
    fitness_log = str(fitness_results.get("fitness_log"))
    diversity_log = str(fitness_results.get("diversity_log"))
    mutation_rate_log = str(fitness_results.get("mutation_rate_log"))
    crossover_rate_log = str(fitness_results.get("crossover_rate_log"))

    # format mutation methods and weights into single string
    #mutation_methods = "".join(f"{method}:{weight}" for method, weight in hyperparameters.get("mutation_weights"))

    csv_row = {
        # Metadata
        "run_id": run_id,
        "timestamp": timestamp,
        "runtime": runtime,

        # Hyperparameters
        "population_size": hyperparameters["population_size"],
        "num_generations": hyperparameters["num_generations"],
        "carrying_capacity": hyperparameters["carrying_capacity"],
        #"mutation_method": mutation_methods, TODO - broken
        "crossover_methods": hyperparameters["crossover_methods"],
        "selection_method": hyperparameters["selection_method"],
        "elitism_weight": hyperparameters["elitism_weight"],

        # Logs
        "fitness_log": fitness_log,
        "diversity_log": diversity_log,
        "mutation_rate_log": mutation_rate_log,
        "crossover_rate_log": crossover_rate_log,

        # Metrics
        "cpu_usage": round(cpu_usage, 3),
        "population_type": population_type,
        "fitness_config": str(fitness_config)
    }

    file_exists = os.path.exists(file_path)
    with open(file_path, 'a', newline="", encoding="utf-8") as file:
        # quoting enabled since many of the data contain commas
        writer = csv.DictWriter(file, fieldnames=csv_row.keys(), quoting=csv.QUOTE_ALL)
        if not file_exists:
            writer.writeheader()
        writer.writerow(csv_row)