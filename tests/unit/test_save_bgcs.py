import os
import pickle
import unittest
import networkx as nx
from unittest.mock import patch, MagicMock
import tempfile
from src.data_pipeline.load_bgcs import load_bgc_jsons
from src.data_pipeline.save_bgcs import save_bgcs
from src.data_pipeline.mol_to_graph import mol_to_graph

class TestSaveBgcs(unittest.TestCase):

    def setUp(self):
        # Set up temporary directory to hold test data
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_file = os.path.join(self.temp_dir.name, "test_preprocessed_bgcs.pkl")

        """
        Mock BGC data, the first compound in the 3rd BGC has been manually
        invalidated for testing error handling functionality.
        """
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
                        "structure": "O=C(N[C@@H]C(N[C@@H](CCN)C(N[C@H](CC2=CC=CC=C2)C(N[C@@H]([C@H](O)C)C(N[C@@H](CCN)C(N[C@H](C(N[C@H]1[C@H](O)C)=O)CCN)=O)=O)=O)=O)=O)[C@@H](CCN)NC([C@@H](NC([C@@H](NC(CCCCC(C)CC)=O)CCN)=O)[C@H](O)C)=O",
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

    # Mock dependencies
    @patch("src.data_pipeline.save_bgcs.load_bgc_jsons")
    @patch("src.data_pipeline.mol_to_graph.mol_to_graph")
    def test_save_bgcs(self, mock_mol_to_graph, mock_load_bgc_jsons):
        # Test that BGCs are properly saved to pickle file

        # Mock load_bgc_json output / BGC data
        mock_load_bgc_jsons.return_value = self.mock_bgc_json

        # Mock graph building functionality
        def mock_graph(smiles_str):
            """
            Manually defined, valid SMILES strings from the test data.
            Invalidated SMILES string omitted (returns None)
            """
            return "GraphObject" if smiles_str in [
                "CCCC(O[C@H]1C[C@](C)(N)[C@H](O)[C@H](C)O1)C(C)C(O)C(CC)\\C=C\\C(O)C(C)C1C\\C=C(C)\\C(O)C(C)C(CC(O)C(C)C(O)CC2CC(O)C(O)C(O)(CC(O[C@@H]3O[C@H](C)[C@@H](O)[C@H](O[C@H]4C[C@@H](N)[C@H](O)[C@@H](C)O4)[C@H]3O[C@@H]3O[C@H](C)[C@@H](O)[C@H](O)[C@H]3O)C(C)CCC(O)CC(O)C\\C=C(CC)\\C(=O)O1)O2)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O",
                "CCC(C)C(C(=O)OC(/C=C/C=C/C=C/C(=O)O)C1(CO1)C)OC(=O)C(C(C)(C)O)O",
                "O=C(N[C@@H](CCNC1=O)C(N[C@@H](CCN)C(N[C@H](CC2=CC=CC=C2)C(N[C@@H]([C@H](O)C)C(N[C@@H](CCN)C(N[C@H](C(N[C@H]1[C@H](O)C)=O)CCN)=O)=O)=O)=O)=O)[C@@H](CCN)NC([C@@H](NC([C@@H](NC(CCCCC(C)C)=O)CCN)=O)[C@H](O)C)=O"
                ] else None
        
        # Dynamic, mock output
        mock_mol_to_graph.side_effect = mock_graph

        save_bgcs(json_dir=None, output_file=self.output_file)

        # Verify creation of pickle file
        self.assertTrue(os.path.exists(self.output_file))
        

        # Load pickle test file
        with open(self.output_file, "rb") as f:
            pkl_data = pickle.load(f)
        
        # Verify all 3 BGCs are there
        self.assertEqual(len(pkl_data), 3)

        # Verify metadata and graph are present for each BGC
        bgc1 = pkl_data[0]
        self.assertEqual(bgc1["accession"], "BGC0000002")
        self.assertIsInstance(bgc1["compounds"][0]["mol_graph"], nx.Graph)

        bgc2 = pkl_data[1]
        self.assertEqual(bgc2["accession"], "BGC0000003")
        self.assertIsInstance(bgc2["compounds"][0]["mol_graph"], nx.Graph)

        bgc3 = pkl_data[2]
        self.assertEqual(bgc3["accession"], "BGC0003170")
        # Verify manually invalidated Mol doesn't produce a graph
        self.assertIsNone(bgc3["compounds"][0]["mol_graph"])
        self.assertIsInstance(bgc3["compounds"][1]["mol_graph"], nx.Graph)

    def tearDown(self):
        self.temp_dir.cleanup()

if __name__ == "__main__":
    unittest.main()