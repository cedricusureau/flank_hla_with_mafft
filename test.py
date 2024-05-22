import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO

#Merge every fasta file ine data/reconstructed_fasta/ into a single file

def merge_fasta_files(input_files, output_file):
    """
    Merge les fichiers FASTA d'entrée en un seul fichier FASTA de sortie.

    Args:
        input_files (list): Une liste de chemins de fichiers FASTA d'entrée.
        output_file (str): Le chemin du fichier FASTA de sortie.
    """
    records = []
    for input_file in input_files:
        records.extend(SeqIO.parse(input_file, "fasta"))

    SeqIO.write(records, output_file, "fasta")
    print(f"Les fichiers FASTA d'entrée ont été fusionnés avec succès en {output_file}")

# Merge file in combined_inputs
combined_input_files = os.listdir("data/merged_fasta")
combined_input_files = [f"data/merged_fasta/{file}" for file in combined_input_files]
merge_fasta_files(combined_input_files,"data/hla_flanked_10kb.fasta")