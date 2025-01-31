import os
import json

def load_bgc_jsons(json_dir):
    """
    Load in BGC JSON files 

    :param json_dir: Directory contain BGC JSON files.
    :returns: List of parsed BGC data.
    """
    bgc_data = []
    for file in os.listdir(json_dir):
        if file.endswith(".json"):
            file_path = os.path.join(json_dir, file)
            with open(file_path, "r") as f:
                data = json.load(f)
                bgc_data.append(data)

    return bgc_data

