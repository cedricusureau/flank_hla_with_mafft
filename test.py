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

merge_fasta_files(["data/reconstructed_fasta/reconstructed_A_output.fasta",
                   "data/reconstructed_fasta/reconstructed_B_output.fasta",
                   "data/reconstructed_fasta/reconstructed_C_output.fasta",
                   "data/reconstructed_fasta/reconstructed_DPB1_output.fasta",
                   "data/reconstructed_fasta/reconstructed_DQB1_output.fasta",
                   "data/reconstructed_fasta/reconstructed_DRB1_output.fasta",
                   "data/reconstructed_fasta/reconstructed_DRB3_output.fasta",
                   "data/reconstructed_fasta/reconstructed_DRB4_output.fasta",
                   "data/reconstructed_fasta/reconstructed_DRB5_output.fasta",
                   "data/reconstructed_fasta/reconstructed_DPA1_output.fasta",
                   "data/reconstructed_fasta/reconstructed_DQA1_output.fasta"],
                  "data/hla_flanked_10kb.fasta")