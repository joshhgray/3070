import os
import unittest
import pytest #type:ignore
import networkx as nx
from src.controller import start_ga
from src.evolutionary_system.utils.config_loader import load_config
from src.evolutionary_system.utils.logger import log_metrics
from unittest.mock import patch
import psutil

class TestController(unittest.TestCase):
    def setUp(self):
        self.file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../config.yaml")
        self.config = load_config()

    def test_config_file_exists(self):
        """ Test that the config file exists """
        self.assertTrue(os.path.isfile(self.file_path))
        
    def test_load_config(self):
        """ Test that the config file is fully loaded"""
        self.assertIsInstance(self.config, dict)
        self.assertIn("crossover_rate", self.config)
        self.assertIn("fitness_weights", self.config)
        self.assertIn("mutation_rate", self.config)
        self.assertIn("num_elite_groups", self.config)
        self.assertIn("num_elite_individuals", self.config)
        self.assertIn("num_generations", self.config)
        self.assertIn("num_threads", self.config)
        self.assertIn("population_size", self.config)
        self.assertIn("selection_method", self.config)

    # TODO - Test controller execution


if __name__ == "__main__":
    unittest.main()