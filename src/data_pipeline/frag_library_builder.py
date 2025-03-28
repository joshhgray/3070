from rdkit import Chem
from rdkit.Chem import BRICS, rdMolDescriptors
import os
import pickle
from collections import Counter
import csv
from tqdm import tqdm

max_atoms = 100 # For efficiency - BRICS can take a while on large mols

# Load BGCs from pickle
cd = os.path.dirname(os.path.abspath(__file__))
population_pool = os.path.join(cd, "../../preprocessed_bgcs.pkl")
with open(population_pool, "rb") as f:
    bgc_data = pickle.load(f)

# Extract and compile All SMILES strings from BGC data
smiles = []
for bgc in bgc_data:
    for compound in bgc.get("compounds", []):
        structure = compound.get("structure")
        if structure:
            smiles.append(structure)

# Convert SMILES to RDKit Mols
mols = [Chem.MolFromSmiles(structure) for structure in smiles if Chem.MolFromSmiles(structure)]

# Decompose Mols with BRICS
frag_library = []
for mol in tqdm(mols):
    if mol is None:
        continue
    # Skip overly large mols
    atoms = rdMolDescriptors.CalcNumHeavyAtoms(mol)
    if atoms > max_atoms:
        continue
    try:
        frags = BRICS.BRICSDecompose(mol, minFragmentSize=2)
        frag_library.extend(frags)
    except Exception as e:
        #print(f"Error performing BRICS Decomposition: {e}")
        continue

# Extract counts for number of occurences per unique SMILES
frag_counts = Counter(frag_library)

# Filter out small unmeaningful and overly large fragmen ts 
frag_counts_filtered = [
    {"fragment_smiles": smiles, "count": count}
    for smiles, count in frag_counts.items()
    if 6 <= len(smiles) <= 60
]

# Sort by most common
frag_counts_filtered_sorted = sorted(frag_counts_filtered, key=lambda frag: frag["count"], reverse=True)

print("Saving...")
# Save fragment library
with open("bgc_fragment_library.csv", "w") as f:
    writer = csv.DictWriter(f, fieldnames=["fragment_smiles", "count"])
    writer.writeheader()
    # Only saving most frequent 1000
    writer.writerows(frag_counts_filtered_sorted[:1000])

print(f"Extracted {len(frag_counts)} fragmentss.")
print(frag_counts.most_common(10))