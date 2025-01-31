import os
import pickle
from src.data_pipeline.load_bgcs import load_bgc_jsons
from src.data_pipeline.mol_to_graph import mol_to_graph

def save_bgcs(json_dir, output_file):
    """
    Loads BGC JSONs, builds each BGCs respective molecular graph
    and stores everything in a pickle file

    :param json_dir: json_dir of unprocessed BGC data in JSON format
    :param output_file: pickle file to store preprocessed BGCs
    """
    bgc_data = load_bgc_jsons(json_dir)

    for bgc in bgc_data:
        for compound in bgc.get("compounds", []):
            smiles_str = compound.get("structure")
            compound["mol_graph"] = mol_to_graph(smiles_str)
    
    with open(output_file, "wb") as f:
        pickle.dump(bgc_data, f)

if __name__ == "__main__":
    cd = os.path.dirname(os.path.abspath(__file__))
    json_dir = os.path.join(cd, '../../mibig_json_4.0')
    output_file = os.path.join(cd, '../../preprocessed_bgcs.pkl')

    save_bgcs(json_dir, output_file)