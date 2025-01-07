import os
import unittest
import csv

from src.evolutionary_system.utils.logger import log_metrics

class TestLogger(unittest.TestCase):

    def setUp(self):
        # Get path to metrics log
        self.cd = os.path.dirname(os.path.abspath(__file__))
        self.file_path = os.path.join(self.cd, "../../src/evolutionary_system/utils/metrics_log.csv")
        # Resolve the file path
        self.file_path = os.path.normpath(self.file_path) 


    def test_log_metrics(self):
        # Mocks/test inputs
        run_id = 'test_10'
        timestamp = "2025-01-01T00:00:00"
        runtime = 30

        hyperparameters = {
            "population_size": 100,
            "num_generations": 20,
            "mutation_rate": 0.1,
            "crossover_rate": 0.1,
            "num_elite_individuals": 1,
            "num_elite_groups": 1,
            "fitness_weights": [1.0, 1.0],
            "num_threads": 1,
        }

        fitness_results = {
            "diversity": [0.94, 0.93, 0.92],
            # TODO - add other fitness results
        }
        validation_results = {
            "docking_score": 0.1,
            # TODO - Add other metrics
        }

        # mock cpu and memory usage
        mem_usage = 2048
        cpu_usage = 90

        log_metrics(
            run_id=run_id,
            timestamp=timestamp,
            runtime=runtime,
            hyperparameters=hyperparameters,
            fitness_results=fitness_results,
            validation_results=validation_results,
            mem_usage=mem_usage,
            cpu_usage=cpu_usage,
        )
        
        # Verify the file exists
        self.assertTrue(os.path.exists(self.file_path))

        with open(self.file_path, "r") as file:
            reader = csv.DictReader(file)
            last_row = list(reader)[-1]

            # Verify each metric logged in metrics_log
            self.assertEqual(last_row["run_id"], run_id)
            self.assertEqual(last_row["timestamp"], timestamp)
            self.assertEqual(float(last_row["runtime"]), runtime)

            # Verify hyperparameters logged
            # TODO - do

            # Verify fitness and validation metrics logged
            fitness_results = last_row["fitness_results"]
            validation_results = last_row["validation_results"]
            self.assertIn("diversity", fitness_results)
            self.assertIn("docking_score", validation_results)

            # Verify performance metrics logged
            self.assertEqual(float(last_row["mem_usage"]), mem_usage)
            self.assertEqual(float(last_row["cpu_usage"]), cpu_usage)

if __name__ == "__main__":
    unittest.main()