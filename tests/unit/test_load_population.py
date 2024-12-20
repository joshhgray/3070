from src.data_pipeline.load_population  import load_initial_population
import unittest
from rdkit import Chem

class TestLoadPopulation(unittest.TestCase):
    def test_load_initial_population(self):
        initial_population = load_initial_population()
        self.assertTrue(all(Chem.MolFromSmiles(smiles_str) is not None for smiles_str in initial_population))
    
if __name__ == "__main__":
    unittest.main()