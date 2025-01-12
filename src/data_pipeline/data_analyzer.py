from Bio import SeqIO
import os

def analyze_data(file_path, fragment_size=1000):
    """
    Analyze the metagenome data file to estimate the size of the population pool

    :param file_path: Where the metagenome data is stored (.gb file)
    :param fragment_size: Size of fragments in population

    :return: Summary of analysis
    """
    total_sequences = 0
    total_length = 0
    fragment_count = 0

    with open(file_path, "r") as file:
        for record in SeqIO.parse(file, "genbank"):
            seq_length = len(record.seq)
            total_sequences += 1
            total_length += seq_length
            fragment_count += seq_length // fragment_size
    
    avg_length = total_length / total_sequences

    stats = {
        "Total Sequences": total_sequences,
        "Total Length (bp)": total_length,
        "Average Sequence Length (bp)": avg_length,
        "Population Size (total fragments (approx))": fragment_count,
    }
    return stats


cd = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(cd, "../../data/metagenome_data.gb")

# TODO - Decide on a size
fragment_size = 1000

summary = analyze_data(file_path, fragment_size)

for key, value in summary.items():
    print(f"{key}: {value}")