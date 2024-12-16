import unittest

from src.evolutionary_system.utils.logger import log_metrics

class TestLogger(unittest.TestCase):
    def test_log_metrics(self):
        # Mocks/test inputs
        mock_time = [100.00, 124.56]
        run_id = 'test_0',
        hyperparameters = {
            "population": 100,
            "n_generations": 20,
            # test additional hyperparameters
        }
        fitness_results = {
            "population_diversity": 10,
            # test additional fitness results
        }
        validation_results = {
            "accuracy": 0.7
            # test additional validation results
        }
        
        log_metrics(
            run_id=run_id,
            hyperparameters=hyperparameters,
            fitness_results=fitness_results,
            validation_results=validation_results,
            start_time=mock_time[0],
            end_time=mock_time[1],
        )

        with open("metrics_log.csv", "r") as file:
            lines = file.readlines()

        self.assertIn("test_0", lines[-1]) # assert run_id is included
        self.assertIn("2024-12", lines[-1]) # assert timestamp is included

if __name__ == "__main__":
    unittest.main()