import os
import pickle
from src.data_pipeline.load_bgcs import load_bgc_jsons
from rdkit import Chem

def save_bgcs(json_dir, output_file):
    """
    Loads BGC JSONs. Filter out BGCs with missing critical data.
    and save to pickle file.

    :param json_dir: json_dir of unprocessed BGC data in JSON format
    :param output_file: pickle file to store preprocessed BGCs
    """
    bgc_data = load_bgc_jsons(json_dir)

    for bgc in bgc_data:
        for compound in bgc.get("compounds", []):
            smiles_str = compound.get("structure")
            if smiles_str:
                compound["mol"] = Chem.MolFromSmiles(smiles_str)
    
    filtered_bgc_data = [
        bgc for bgc in bgc_data
        if any(c.get("mol") is not None for c in bgc.get("compounds", []))
    ]

    with open(output_file, "wb") as f:
        pickle.dump(filtered_bgc_data, f)
    
    print(f"Successfully loaded {len(filtered_bgc_data)} BGCs to the population pool.")

if __name__ == "__main__":
    cd = os.path.dirname(os.path.abspath(__file__))
    json_dir = os.path.join(cd, '../../mibig_json_4.0')
    output_file = os.path.join(cd, '../../preprocessed_bgcs.pkl')

    save_bgcs(json_dir, output_file)