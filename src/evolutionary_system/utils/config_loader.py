import os
import yaml

def load_config():
    file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../config.yaml")
    
    with open(file_path, "r") as file:
        config = yaml.safe_load(file)
        
    return config
