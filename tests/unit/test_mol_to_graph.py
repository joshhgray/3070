import unittest
import networkx as nx
from rdkit import Chem
from src.data_pipeline.mol_to_graph import mol_to_graph

class TestMolToGraph(unittest.TestCase):

    def setUp(self):
        self.valid_smiles_strs = [
            # Random smiles strings sourced from: https://www.wolframcloud.com/env/a7fc205c-066b-45c2-bea6-90069bd09da4?src=CloudBasicCopiedContent#sidebar=basic-notebook-links
            "c1([N+](=O)[O-])c[nH]c(Cl)c1Cl",
            "[nH]1cc(c(Cl)c1Cl)[N+](=O)[O-]",
            "c1[nH]c(Cl)c(c1[N+]([O-])=O)Cl",
            "c1([nH]cc(c1Cl)[N+]([O-])=O)Cl",
            "Clc1c(c(c[nH]1)[N+]([O-])=O)Cl"
        ]
        self.invalid_smiles_strs = [
            # Manually invalidated
            "assdddf", # bond failure
            "sxxfds", # dangling ring closure
            "",
            None,
            0
        ]
    def test_valid_smiles_str(self):
        for smiles_str in self.valid_smiles_strs:
            graph = mol_to_graph(smiles_str)
            self.assertIsInstance(graph, nx.Graph)
            self.assertGreater(len(graph.nodes), 0)
            self.assertGreaterEqual(len(graph.edges), 0)

    def test_invalid_smiles_str(self):
        for invalid_smiles_str in self.invalid_smiles_strs:
            graph = mol_to_graph(invalid_smiles_str)
            print(graph)
            self.assertIsNone(graph)
    
if __name__ == "__main__":
    unittest.main()