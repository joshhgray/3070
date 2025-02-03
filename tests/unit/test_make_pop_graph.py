import unittest
import networkx as nx
from unittest.mock import patch
from src.data_pipeline.make_population_graph import make_population_graph

class TestMakePopGraph(unittest.TestCase):
    def setUp(self):
        # Mock BGC dataset
        self.mock_bgc_json = [
            {
                "accession": "BGC0000002",
                "biosynthesis": {
                    "classes": [
                        {
                            "class": "PKS",
                        }
                    ]
                },
                "compounds": [
                    {
                        "name": "aculeximycin",
                        "structure": "CCCC(O[C@H]1C[C@](C)(N)[C@H](O)[C@H](C)O1)C(C)C(O)C(CC)\\C=C\\C(O)C(C)C1C\\C=C(C)\\C(O)C(C)C(CC(O)C(C)C(O)CC2CC(O)C(O)C(O)(CC(O[C@@H]3O[C@H](C)[C@@H](O)[C@H](O[C@H]4C[C@@H](N)[C@H](O)[C@@H](C)O4)[C@H]3O[C@@H]3O[C@H](C)[C@@H](O)[C@H](O)[C@H]3O)C(C)CCC(O)CC(O)C\\C=C(CC)\\C(=O)O1)O2)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O",
                        "mass": 1672.9651350679992,
                        "formula": "C81H144N2O33",
                        "mol_graph": nx.Graph()
                    }
                ],
                "taxonomy": {
                    "name": "Kutzneria albida DSM 43870",
                },
            },
            {
                "accession": "BGC0000003",
                "biosynthesis": {
                    "classes": [
                        {
                            "class": "PKS",
                        }
                    ]
                },
                "compounds": [
                    {
                        "structure": "CCC(C)C(C(=O)OC(/C=C/C=C/C=C/C(=O)O)C1(CO1)C)OC(=O)C(C(C)(C)O)O",
                        "mass": 440.20463260399987,
                        "formula": "C22H32O9",
                        "mol_graph": nx.Graph()
                    }
                ],
                "taxonomy": {
                    "name": "Alternaria alternata",
                },
            },
            {
                "accession": "BGC0003170",
                "biosynthesis": {
                    "classes": [
                        {
                            "class": "NRPS",
                        }
                    ]
                },
                "compounds": [
                    {
                        "structure": "O=C(N[C@@H](CCNC1=O)C(N[C@@H](CCN)C(N[C@H](CC2=CC=CC=C2)C(N[C@@H]([C@H](O)C)C(N[C@@H](CCN)C(N[C@H](C(N[C@H]1[C@H](O)C)=O)CCN)=O)=O)=O)=O)=O)[C@@H](CCN)NC([C@@H](NC([C@@H](NC(CCCCC(C)CC)=O)CCN)=O)[C@H](O)C)=O",
                        "mass": 1190.7135,
                        "formula": "C54H94N16O14",
                        "mol_graph": nx.Graph()
                    },
                    {
                        "structure": "O=C(N[C@@H](CCNC1=O)C(N[C@@H](CCN)C(N[C@H](CC2=CC=CC=C2)C(N[C@@H]([C@H](O)C)C(N[C@@H](CCN)C(N[C@H](C(N[C@H]1[C@H](O)C)=O)CCN)=O)=O)=O)=O)=O)[C@@H](CCN)NC([C@@H](NC([C@@H](NC(CCCCC(C)C)=O)CCN)=O)[C@H](O)C)=O",
                        "mass": 1176.6979,
                        "formula": "C53H92N16O14",
                        "mol_graph": None
                    }
                ],
                "taxonomy": {
                    "name": "Paenibacillus polymyxa M1",
                    "ncbiTaxId": 1052684
                }
            }
        ]
    
    @patch("src.data_pipeline.make_population_graph.sample_bgcs")
    def test_make_population_graph(self, mock_sample_bgcs):
        mock_sample_bgcs.return_value = self.mock_bgc_json

        # Make a 3-BGC graph for testings
        test_population_graph = make_population_graph(3)

        # Verify root node exists
        self.assertIn("Population", test_population_graph.nodes)
        self.assertEqual(test_population_graph.nodes["Population"]["level"], "Root")


        # Verify group nodes exists
        group_nodes = [node for node in test_population_graph.nodes if test_population_graph.nodes[node]["level"] == "Group"]

        self.assertEqual(len(group_nodes), 2)
        self.assertIn("PKS", group_nodes)
        self.assertIn("NRPS", group_nodes)

        # Verify inidivudal nodes exist
        individual_nodes = [node for node in test_population_graph.nodes if test_population_graph.nodes[node]["level"] == "Individual"]
        self.assertEqual(len(individual_nodes), 3)


        # Verify Population to Group Edges
        population_edges = list(test_population_graph.edges("Population"))

        # Verify proper splitting of groups
        self.assertEqual(len(population_edges), 2)

        # Verify Population to Group Edges
        self.assertIn(("Population", "PKS"), population_edges)
        self.assertIn(("Population", "NRPS"), population_edges)

        # Verify Group to Individual Edges
        pks_edges = list(test_population_graph.edges("PKS"))
        nrps_edges = list(test_population_graph.edges("NRPS"))

        # Verify edges exist
        self.assertEqual(len(pks_edges), 2)
        self.assertEqual(len(nrps_edges), 1)

        # Validate connections
        self.assertIn(("PKS", "BGC0000002"), pks_edges)
        self.assertIn(("PKS", "BGC0000003"), pks_edges)
        self.assertIn(("NRPS", "BGC0003170"), nrps_edges)

        # Verify that the two compounds are present
        multi_compound_bgc = test_population_graph.nodes["BGC0003170"]

        self.assertEqual(len(multi_compound_bgc["compounds"]), 2)
        
        # Verify contents of compounds
        compound1 = multi_compound_bgc["compounds"][0]
        compound2 = multi_compound_bgc["compounds"][1]

        self.assertEqual(compound1["mass"], 1190.7135)
        self.assertEqual(compound2["mass"], 1176.6979)
        self.assertEqual(compound1["formula"], "C54H94N16O14")
        self.assertEqual(compound2["formula"], "C53H92N16O14")

        # Verify mol_graphs are handled properly
        bgc1 = test_population_graph.nodes["BGC0000002"]
        bgc2 = test_population_graph.nodes["BGC0000003"]
        bgc3 = test_population_graph.nodes["BGC0003170"]

        self.assertIsInstance(bgc1["compounds"][0]["mol_graph"], nx.Graph)
        self.assertIsInstance(bgc2["compounds"][0]["mol_graph"], nx.Graph)
        self.assertIsInstance(bgc3["compounds"][0]["mol_graph"], nx.Graph)
        self.assertIsNone(bgc3["compounds"][1]["mol_graph"])


if __name__ == "__main__":
    unittest.main()