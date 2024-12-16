import unittest
from src.data_pipeline.make_graph import make_graph

class TestMakeGraph(unittest.TestCase):
    def setUp(self):
        # Random SMILES strings
        self.test_population = [
            "CC(CO)=CCNC3(=NC=NC2(N(C1(C(C(C(O1)COP([O-])([O-])=O)O)O))C=NC=23))",
            "CCOC(=O)C(C(=O)OCC)N.Cl",
            "CC5(C=CC4(C(=CC=C1(C(OCC2(O)(C3(C=CC(O)=CC(O[CH]12)=3)))=4))O5))(C)",
            "CCOP(=O)(CC(=O)O)OCC", 
            "C2(CC1(C=C(OC)C(=CC=1[CH]3(CC4(=CC=C(OC)C(OC)=C(CN23)4)))OC))",
            "COC(=O)COC1=C(C=C(C=C1)Cl)Cl",
            "C1=COC(=C1)CN",
            "COC(=O)C1=C(C=CC(=C1)[N+](=O)[O-])Cl",
            "C(=C(F)OC(C(C(F)(F)F)(F)F)(F)F)(F)F.C(=C(F)F)(F)F",
            "NC1(CC1)P(O)([O-])=O",
            "C1=CC=C(C(=C1)C(=S)N)Cl",
            "CCCC1=NC=CN=C1C",
            "COC(=O)C1=CC(=O)NO1",
            "COC1=CC=C(C=C1)C2=CC=C(C=C2)Br",
            "COC1(=CC(=C(C=C1)C2(C(C3(=C(C=C(C=C(OC2)3)OC4(C(C(C(C(O4)CO)O)O)O))O))=O))OC)",
            "CC1=CC=C(C=C1)S(=O)(=O)N(C)N=O",
            "CC(=O)OC(C1=CC=C(O1)[N+](=O)[O-])OC(=O)C",
            "CCCC[N+](C)(CCCC)CCCC.CCCCOP(=O)([O-])OCCCC",
            "C[SiH](C)O[Si](C)(C)O[SiH](C)C",
            "CCOC(=O)C1CCCNC1",
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCOC(=O)CCCCCCCCCCCCCCCCCCCCC",
            "CCCCCCCC(=O)[O-].[Na+]",
            "CC1(C)(CC5(C(CC1)(CCC4(C)(C3(C)(C(C2(C)(C(C(C)(C)C(O)C(O)C2)CC3))CC=C45)))C([O-])=O))",
            "CC(=CC([O-])=O)C=CC1(O)(C(C)=CC(=O)CC(CO)(C)1)",
            "CC1=C(NC(SC1)C(NC(=O)CC1=CC=CS1)C(O)=O)C(O)=O",
            "C(CSSCCN)N.Cl.Cl",
            "NC1=C2N=CN(C3OC(COP(O)(=O)OP(O)(=O)OCC4OC(O)C(O)C4O)C(O)C3O)C2=NC=N1",
            "CC(=O)OC(C)(C)CC1=CC=CC=C1",
            "CC(C)(C)C(=O)NC1=CC=CC=C1Cl",
            "C=CC(=O)OCC.C=C",
            "[Cl-].[Cl-].[Dy+2]",
            "C1C(C(OC1N2C=CC(=NC2=O)N)COP(=O)(O)O)O",
            "COC(CO)C1=CC=CC=C1",
            "CC1=C(C(=C(C1[Si](C)(C)NC(C)(C)C)C)C)C",
            "COC(=O)CN1C=CN=C1[N+](=O)[O-]",
            "CCCCCCCCC1CC[NH+](C)C1.C(F)(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F",
            "CC(C)C1=NC(=NC(=C1C=CC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)N(C)S(=O)(=O)C",
            "C1=CC=NC=C1.C(F)(F)(F)S(=O)(=O)O",
            "CC(=O)CC(C1=CC=C(C=C1)[N+](=O)[O-])C2=C(OC3=CC=CC=C3C2=O)O",
            "C1COC(O1)CC[P+](C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=CC=C4.[Br-]",
            "CC(=O)NC(C(=O)O)C(C)(C)S","C1=NC(=O)NC2=C1N(C=N2)C3C(C(C(O3)COP(=O)(O)O)O)O",
            "CC1=C(C(C(=C(N1)C)C(=O)OCC(C)C)C2=CC=CC=C2[N+](=O)[O-])C(=O)OC",
            "CCC#CC1=CC=CC=C1",
            "CC(C)(C)C(CO)NC=C1C=C(C=C(C1=O)I)I",
            "C(C(=O)O)OP(=O)(O)O"
        ]
    
    def test_graph_exists(self):
        # Test that the components of the graph exist
        
        graph = make_graph(self.test_population)

        # Test root node
        self.assertIn("Population", graph.nodes)
        # Test group nodes

        # test individual nodes

if __name__ == "__main__":
    unittest.main()
