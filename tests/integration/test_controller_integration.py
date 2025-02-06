import os
import unittest
import pytest #type:ignore
import networkx as nx
from src.controller import start_ga
from src.evolutionary_system.utils.config_loader import load_config

class TestControllerIntegration(unittest.TestCase):
    def setUp(self):
        self.file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../config.yaml")
        self.config = load_config()


    def test_controller_to_ga(self):
        # Test controller integration with GA
        final_population, diversity_log = start_ga()

        # Verify final population
        self.assertIsInstance(final_population, nx.DiGraph)

        # Verify diversity log
        self.assertEqual(len(diversity_log), self.config["num_generations"])

if __name__ == "__main__":
    unittest.main()