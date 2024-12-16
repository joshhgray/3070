import unittest

from src.evolutionary_system.fitness import calculate_qed

class TestFitness(unittest.TestCase):
    def setUp(self):
        # Random SMILES string
        self.smiles_str = "CC(CO)=CCNC3(=NC=NC2(N(C1(C(C(C(O1)COP([O-])([O-])=O)O)O))C=NC=23))"

    def testQEDFunction(self):
        # Test that the fitness function returns a valid score

        fitness_score = calculate_qed(self.smiles_str)
        self.assertNotEqual(fitness_score, 0)

if __name__ == "__main__":
    unittest.main()