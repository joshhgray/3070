import csv
import os

def log_metrics(run_id, time hyperparameters, fitness_results, validation_results):
    # connect to metrics log
    cd = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(cd, "metrics_log.csv")
    
    


    # get hyperparameters (e.g., population size, n_generations)
    # get fitness/validation results

    csv_row = {
        "timestamp": timestamp,
        "run_id": run_id,
        "runtime": runtime,
        # write out hyperparameters
        # write out fitness/validation results
        # write out performance metrics
    }
    file_exists = os.path.exists(file_path)
    with open(file_path, 'a', newline="") as file:
        writer = csv.DictWriter(file, fieldnames=csv_row.keys())
        if not file_exists:
            writer.writeheader()
        writer.writerow(csv_row)