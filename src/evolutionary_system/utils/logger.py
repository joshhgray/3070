import datetime
import psutil # type: ignore
import csv
import os

def log_metrics(run_id, hyperparameters, fitness_results, validation_results, start_time, end_time):
    # connect to metrics log
    cd = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(cd, "metrics_log.csv")
    
    timestamp = datetime.datetime.now().isoformat()
    runtime = end_time - start_time

    process = psutil.Process()
    # Calcuate starting Memory
    # Calculate starting CPU usage

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