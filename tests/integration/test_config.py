import os
import unittest
from src.evolutionary_system.utils.config_loader import load_config

class TestConfigIntegration(unittest.TestCase):
    def setUp(self):
        self.file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../config.yaml")

        
    def test_config_file_exists(self):
        self.assertTrue(os.path.isfile(self.file_path))
        
    def test_load_config(self):
        config = load_config()
        self.assertIsInstance(config, dict)
        
if __name__ == "__main__":
    unittest.main()