import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO

alignment = AlignIO.read("data/aligned_outputs/cleaned_A_output.clustal", "clustal")

breakpoint()