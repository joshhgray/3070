import unittest
from src.evolutionary_system.fitness import calculate_qed

class TestFitness(unittest.TestCase):
    def setUp(self):
        # Random SMILES string
        self.structure = "CC(CO)=CCNC3(=NC=NC2(N(C1(C(C(C(O1)COP([O-])([O-])=O)O)O))C=NC=23))"

    def test_calculate_qed(self):
        # Test that the fitness function returns a valid score
        self.fitness_score = calculate_qed(self.structure)
        self.assertNotEqual(self.fitness_score, 0)
    


    # TODO - Test other fitness function calculations as they are added
        

if __name__ == "__main__":
    unittest.main()