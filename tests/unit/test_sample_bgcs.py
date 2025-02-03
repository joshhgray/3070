import os
import pickle
import unittest
import tempfile
from unittest.mock import patch, MagicMock
from src.data_pipeline.sample_population import sample_bgcs

class TestSampleBgcs(unittest.TestCase):

    def setUp(self):
        # temp pickle file with mock BGCs for testing
        self.temp_dir = tempfile.TemporaryDirectory()
        self.pickle_file = os.path.join(self.temp_dir.name, 'test_bgcs.pkl')

        self.mock_bgcs = [
            {"accession": "BGC0001231", "biosynthesis": {"classes": [{"class": "NRPS"}]}, "compounds" : {"structure":"O=C([C@@H](O)[C@H](NC(C[C@@H](O)CNC(CNC([C@@H](CC1=CNC2=C1C=CC(Cl)=C2)NC(CN(C(C3)=O)C)=O)=O)=O)=O)[C@H](O)[C@@H](O)C/C=C/C4=CC=CC=C4)N[C@H](C5)[C@@]3(O)NC5=O", "mass": 868.3158466880002} },
            {"accession": "BGC0003170", "biosynthesis": {"classes": [{"class": "NRPS"}]}, "compounds" : {"structure":"O=C(N[C@@H](CCNC1=O)C(N[C@@H](CCN)C(N[C@H](CC2=CC=CC=C2)C(N[C@@H]([C@H](O)C)C(N[C@@H](CCN)C(N[C@H](C(N[C@H]1[C@H](O)C)=O)CCN)=O)=O)=O)=O)=O)[C@@H](CCN)NC([C@@H](NC([C@@H](NC(CCCCC(C)CC)=O)CCN)=O)[C@H](O)C)=O", "mass": 1190.7135} },
            {"accession": "BGC0001388", "biosynthesis": {"classes": [{"class": "ribosomal"}]}, "compounds" : {"structure":"invalid", "mass": 0} }
        ]

        with open(self.pickle_file, "wb") as f:
            pickle.dump(self.mock_bgcs, f)

    def test_sample_bgcs(self):
        test_pop = sample_bgcs(self.pickle_file, population_size=2)
        self.assertEqual(len(test_pop), 2)
        self.assertTrue(all(bgc in self.mock_bgcs for bgc in test_pop))

    def test_sample_bgcs_invalid_input(self):
        test_pop = sample_bgcs(self.pickle_file, population_size=40)
        self.assertEqual(len(test_pop), 3)

if __name__ == "__main__":
    unittest.main()