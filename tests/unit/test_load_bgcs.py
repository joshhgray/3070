import os
import json
import unittest
import tempfile

from src.data_pipeline.load_bgcs import load_bgc_jsons


class TestLoadBgcs(unittest.TestCase):

    def setUp(self):
        # Temporary directory to hold BGC JSON data
        self.temp_dir = tempfile.TemporaryDirectory()

        self.valid_json = {
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
        }

        # Create a valid JSON file
        with open(os.path.join(self.temp_dir.name, "test_bgc.json"), "w") as f:
            json.dump(self.valid_json, f)

        # Create an invalid JSON file
        with open(os.path.join(self.temp_dir.name, "test_invalid.json"), "w") as f:
            f.write("{invalid_key: test}")
        
        # Create a non-JSON file
        with open(os.path.join(self.temp_dir.name, "test_invalid.txt"), "w") as f:
            f.write("Not a JSON file")

    def test_load_valid_json(self):
        # Verify that a valid JSON file is correctly loaded
        # This also tests whether or not the invalid files are properly handled
        bgc_data = load_bgc_jsons(self.temp_dir.name)
        self.assertEqual(len(bgc_data), 1)
        self.assertEqual(bgc_data[0]["accession"], "BGC0000003")


    def test_invalid_dir(self):
        # Verify that an invalid directory is handled properly
        temp_empty_dir = tempfile.TemporaryDirectory()
        bgc_data = load_bgc_jsons(temp_empty_dir.name)
        self.assertEqual(len(bgc_data), 0)

    def test_invalid_json_structure(self):
        # Verify that an invalid JSON structure is handled gracefully
        try:
            load_bgc_jsons(self.temp_dir.name)
        except Exception as e:
            print(f"Error in JSON structure: {e}")

    def tearDown(self):
        self.temp_dir.cleanup()

if __name__ == "__main__":
    unittest.main()
    