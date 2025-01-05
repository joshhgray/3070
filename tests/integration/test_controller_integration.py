import os
import unittest
import pytest #type:ignore
from src.controller import start_ga
from src.evolutionary_system.utils.config_loader import load_config

class TestControllerIntegration(unittest.TestCase):
    def setUp(self):
        self.file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../config.yaml")

    def test_config_file_exists(self):
        """ Test that the config file exists """
        self.assertTrue(os.path.isfile(self.file_path))
        
    def test_load_config(self):
        """ Test that the config file is fully loaded"""
        config = load_config()
        self.assertIsInstance(config, dict)
        self.assertIn("crossover_rate", config)
        self.assertIn("fitness_weights", config)
        self.assertIn("mutation_rate", config)
        self.assertIn("num_elite_groups", config)
        self.assertIn("num_elite_individuals", config)
        self.assertIn("num_generations", config)
        self.assertIn("num_threads", config)
        self.assertIn("population_size", config)
        self.assertIn("selection_method", config)


    def test_build_initial_population():
        """ 
        Test the initial population is properly loaded 
        and made into a graph.
        """
        # TODO - test whether initial population is made
        assert True # TODO - temp

    def test_controller_to_ga():
        """ Test controller integration with ga """
        final_population, diversity_log = start_ga()

        # TODO - validate final population
        # TODO - validate diversity log
        assert True # TODO - temp

    def test_controller_to_dash():
        """ Test Controller interaction with Dash """
        # TODO -add debug logs
        assert True # TODO - temp
    
if __name__ == "__main__":
    unittest.main()