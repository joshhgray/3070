from Bio import SeqIO #type:ignore
import os

def split_sequences(input_file, output_dir):
    """
    Splits up a single file of multiple gene sequences into individual 
    gene sequence files

    :param input_file: single file containing multilple gene sequcnes
    :param output_dir: location to save individual gene sequence files
    """

    with open(input_file, "r") as file:
        for sequence in SeqIO.parse(file, "genbank"):
            output_file = os.path.join(output_dir, f"{sequence.id}.gb")
            SeqIO.write(sequence, output_file, "genbank")
    

cd = os.path.dirname(os.path.abspath(__file__))

input_file = os.path.join(cd, "../../data/metagenome_data.gb")
output_dir = os.path.join(cd, "../../data/split_sequences")
split_sequences(input_file, output_dir)
