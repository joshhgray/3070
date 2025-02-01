import unittest
from src.data_pipeline.make_population_graph import make_population_graph


class TestMakePopGraph(unittest.TestCase):
    def setUp(self):
        # Mock BGC 
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
                        "formula": "C81H144N2O33"
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
                        "formula": "C22H32O9"
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
                        "formula": "C54H94N16O14"
                    },
                    {
                        "structure": "O=C(N[C@@H](CCNC1=O)C(N[C@@H](CCN)C(N[C@H](CC2=CC=CC=C2)C(N[C@@H]([C@H](O)C)C(N[C@@H](CCN)C(N[C@H](C(N[C@H]1[C@H](O)C)=O)CCN)=O)=O)=O)=O)=O)[C@@H](CCN)NC([C@@H](NC([C@@H](NC(CCCCC(C)C)=O)CCN)=O)[C@H](O)C)=O",
                        "mass": 1176.6979,
                        "formula": "C53H92N16O14"
                    }
                ],
                "taxonomy": {
                    "name": "Paenibacillus polymyxa M1",
                    "ncbiTaxId": 1052684
                }
        }
        ]
    
    def test_graph_nodes(self):
        test_population_graph = make_population_graph(self.mock_bgc_json)

        # Verify root node exists
        self.assertIn("Population", test_population_graph.nodes)
        self.assertEqual(test_population_graph.nodes["Population"]["level"], "Root")

        group_nodes = [node for node in test_population_graph.nodes if test_population_graph.nodes[node]["level"] == "Group"]
        
        # Verify group nodes exists
        self.assertEqual(len(group_nodes), 2)
        self.assertIn("PKS", group_nodes)
        self.assertIn("NRPS", group_nodes)

        individual_nodes = [node for node in test_population_graph.nodes if test_population_graph.nodes[node]["level"] == "Individual"]
        
        # Verify inidivudal nodes exist
        self.assertEqual(len(individual_nodes), 3)

    def test_graph_edges(self):
        test_population_graph = make_population_graph(self.mock_bgc_json)

        # Population to Group Edges
        population_edges = list(test_population_graph.edges("Population"))
        # Should split into 2 groups
        self.assertEqual(len(population_edges), 2)
        # Verify Population to Group Edges
        self.assertIn(("Population", "PKS"), population_edges)
        self.assertIn(("Population", "NRPS"), population_edges)

        # Group to Individual Edges
        pks_edges = list(test_population_graph.edges("PKS"))
        nrps_edges = list(test_population_graph.edges("NRPS"))
        # Verify edges exist
        self.assertEqual(len(pks_edges), 2)
        self.assertEqual(len(nrps_edges), 1)
        # Validate connections
        self.assertIn(("PKS", "BGC0000002"), pks_edges)
        self.assertIn(("PKS", "BGC0000003"), pks_edges)
        self.assertIn(("NRPS", "BGC0003170"), nrps_edges)

        


if __name__ == "__main__":
    unittest.main()