import unittest
import networkx as nx
import rdkit
from rdkit import Chem
from src.evolutionary_system.fitness_operations.fitness_functions import calculate_qed

class TestCalculateQED(unittest.TestCase):
    def setUp(self):
        # Build mock graph object
        self.compound = nx.Graph()

        # Test Compound data
        ATOMS = [(0, {'element': 'C', 'atomic_num': 6, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(4), 'num_explicit_hs': 0, 'is_aromatic': False}), (1, {'element': 'C', 'atomic_num': 6, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(4), 'num_explicit_hs': 0, 'is_aromatic': False}), (2, {'element': 'C', 'atomic_num': 6, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(4), 'num_explicit_hs': 0, 'is_aromatic': False}), (3, {'element': 'C', 'atomic_num': 6, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(4), 'num_explicit_hs': 0, 'is_aromatic': False}), (4, {'element': 'C', 'atomic_num': 6, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(4), 'num_explicit_hs': 0, 'is_aromatic': False}), (5, {'element': 'C', 'atomic_num': 6, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(4), 'num_explicit_hs': 0, 'is_aromatic': False}), (6, {'element': 'O', 'atomic_num': 8, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(3), 'num_explicit_hs': 0, 'is_aromatic': False}), (7, {'element': 'C', 'atomic_num': 6, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(3), 'num_explicit_hs': 0, 'is_aromatic': False}), (8, {'element': 'O', 'atomic_num': 8, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(3), 'num_explicit_hs': 0, 'is_aromatic': False}), (9, {'element': 'C', 'atomic_num': 6, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(3), 'num_explicit_hs': 0, 'is_aromatic': False}), (10, {'element': 'C', 'atomic_num': 6, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(3), 'num_explicit_hs': 0, 'is_aromatic': False}), (11, {'element': 'C', 'atomic_num': 6, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(3), 'num_explicit_hs': 0, 'is_aromatic': False}), (12, {'element': 'C', 'atomic_num': 6, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(3), 'num_explicit_hs': 0, 'is_aromatic': False}), (13, {'element': 'O', 'atomic_num': 8, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(3), 'num_explicit_hs': 0, 'is_aromatic': False}), (14, {'element': 'C', 'atomic_num': 6, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(4), 'num_explicit_hs': 0, 'is_aromatic': False}), (15, {'element': 'C', 'atomic_num': 6, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(4), 'num_explicit_hs': 0, 'is_aromatic': False}), (16, {'element': 'C', 'atomic_num': 6, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(4), 'num_explicit_hs': 0, 'is_aromatic': False}), (17, {'element': 'C', 'atomic_num': 6, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(1), 'hybridization': rdkit.Chem.rdchem.HybridizationType(4), 'num_explicit_hs': 1, 'is_aromatic': False}), (18, {'element': 'O', 'atomic_num': 8, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(4), 'num_explicit_hs': 0, 'is_aromatic': False}), (19, {'element': 'C', 'atomic_num': 6, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(4), 'num_explicit_hs': 0, 'is_aromatic': False}), (20, {'element': 'C', 'atomic_num': 6, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(3), 'num_explicit_hs': 0, 'is_aromatic': False}), (21, {'element': 'O', 'atomic_num': 8, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(3), 'num_explicit_hs': 0, 'is_aromatic': False}), (22, {'element': 'C', 'atomic_num': 6, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(3), 'num_explicit_hs': 0, 'is_aromatic': False}), (23, {'element': 'O', 'atomic_num': 8, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(3), 'num_explicit_hs': 0, 'is_aromatic': False}), (24, {'element': 'C', 'atomic_num': 6, 'formal_charge': 0, 'chiral_tag': rdkit.Chem.rdchem.ChiralType(0), 'hybridization': rdkit.Chem.rdchem.HybridizationType(3), 'num_explicit_hs': 0, 'is_aromatic': False})]
        BONDS = [(0, 1, {'bond_type': 'SINGLE'}), (1, 2, {'bond_type': 'SINGLE'}), (2, 3, {'bond_type': 'SINGLE'}), (3, 4, {'bond_type': 'SINGLE'}), (4, 5, {'bond_type': 'SINGLE'}), (5, 6, {'bond_type': 'SINGLE'}), (5, 24, {'bond_type': 'SINGLE'}), (6, 7, {'bond_type': 'SINGLE'}), (7, 8, {'bond_type': 'DOUBLE'}), (7, 9, {'bond_type': 'SINGLE'}), (9, 10, {'bond_type': 'DOUBLE'}), (10, 11, {'bond_type': 'SINGLE'}), (10, 24, {'bond_type': 'SINGLE'}), (11, 12, {'bond_type': 'DOUBLE'}), (11, 19, {'bond_type': 'SINGLE'}), (12, 13, {'bond_type': 'SINGLE'}), (12, 20, {'bond_type': 'SINGLE'}), (13, 14, {'bond_type': 'SINGLE'}), (14, 15, {'bond_type': 'SINGLE'}), (14, 16, {'bond_type': 'SINGLE'}), (14, 17, {'bond_type': 'SINGLE'}), (17, 18, {'bond_type': 'SINGLE'}), (17, 19, {'bond_type': 'SINGLE'}), (20, 21, {'bond_type': 'DOUBLE'}), (20, 22, {'bond_type': 'SINGLE'}), (22, 23, {'bond_type': 'SINGLE'}), (22, 24, {'bond_type': 'DOUBLE'})]
        
        # Extract and add Atoms and their attributes
        for atom_idx, attributes in ATOMS:
            self.compound.add_node(atom_idx,
                                   element=attributes['element'],
                                   atomic_num=attributes['atomic_num'],
                                   formal_charge=attributes['formal_charge'],
                                   chiral_tag=attributes['chiral_tag'],                                   
                                   hybridization=attributes['hybridization'],
                                   num_explicit_hs=attributes['num_explicit_hs'],
                                   is_aromatic=attributes['is_aromatic']
                                   )
        # Extract and add Bonds and its type
        for atom1, atom2, bond_attributes in BONDS:
            self.compound.add_edge(atom1, atom2, bond_type=bond_attributes['bond_type'])

    def test_calculate_qed(self):

        #print(f"atoms: {self.compound.nodes(data=True)}")
        #print(f"bonds: {self.compound.edges(data=True)}")


        # Verify that the QED calculation returns a valid score
        self.qed_score = calculate_qed(self.compound)
        print(f"QED Score: {self.qed_score}")
        self.assertGreaterEqual(self.qed_score, 0)
        self.assertLessEqual(self.qed_score, 1)
        

if __name__ == "__main__":
    unittest.main()