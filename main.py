import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO

def run_mafft(powershell_script, input_file, output_file, auto=True, clustalout=True, inputorder=True):
    """
    Exécute MAFFT en utilisant un script PowerShell.

    Args:
        powershell_script (str): Le chemin absolu du script PowerShell mafft-signed.ps1.
        input_file (str): Le chemin absolu du fichier d'entrée.
        output_file (str): Le chemin absolu du fichier de sortie.
        auto (bool): Utiliser l'option --auto pour MAFFT.
        clustalout (bool): Utiliser l'option --clustalout pour MAFFT.
        inputorder (bool): Utiliser l'option --inputorder pour MAFFT.
    """
    command = [
        "powershell",
        "-ExecutionPolicy", "Bypass",
        "-File", powershell_script
    ]

    if auto:
        command.append("--auto")
    if clustalout:
        command.append("--clustalout")
    if inputorder:
        command.append("--inputorder")

    command.append(input_file)

    with open(output_file, 'w') as outfile:
        try:
            result = subprocess.run(command, check=True, text=True, stdout=outfile, stderr=subprocess.PIPE)
            print(f"MAFFT exécuté avec succès. La sortie est enregistrée dans {output_file}")
        except subprocess.CalledProcessError as e:
            print(f"Erreur lors de l'exécution de MAFFT: {e.stderr}")

def prepare_and_run_mafft(parameters, powershell_script):
    for locus, (transcript_id, orientation) in parameters.items():
        ensembl_file = f"data/reference_transcripts/HLA-{locus}_10kb.fasta"
        imgt_file = f"data/hla_to_flank/filtered_{locus}_gen_30_mismatch.fasta"
        combined_input_file = f"data/combined_inputs/combined_{locus}_input.fasta"
        clustal_output_file = f"data/aligned_outputs/aligned_{locus}_output.clustal"
        cleaned_clustal_file = f"data/aligned_outputs/cleaned_{locus}_output.clustal"
        fasta_output_file = f"data/reconstructed_fasta/reconstructed_{locus}_output.fasta"
        merged_fasta_file = f"data/merged_fasta/merged_{locus}_output.fasta"

        # Lire les séquences de référence
        ensembl_sequences = SeqIO.to_dict(SeqIO.parse(ensembl_file, "fasta"))
        if transcript_id not in ensembl_sequences:
            print(f"Transcript {transcript_id} not found in {ensembl_file}")
            continue

        reference_seq = ensembl_sequences[transcript_id]
        if orientation == "reverse":
            reference_seq.seq = reference_seq.seq.reverse_complement()

        # Lire les séquences supplémentaires et stocker les en-têtes complets
        imgt_sequences = list(SeqIO.parse(imgt_file, "fasta"))
        imgt_headers = {record.id: record.description for record in imgt_sequences}

        # Ajouter la séquence de référence en premier
        combined_sequences = [reference_seq] + imgt_sequences

        # Écrire les séquences combinées dans un fichier d'entrée
        SeqIO.write(combined_sequences, combined_input_file, "fasta")

        # Appeler MAFFT
        run_mafft(powershell_script, combined_input_file, clustal_output_file)

        # Supprimer les quatre premières lignes du fichier Clustal
        remove_header_lines(clustal_output_file, cleaned_clustal_file)

        # Aligner les séquences au format Clustal
        align_to_reference(cleaned_clustal_file, fasta_output_file, transcript_id, imgt_headers)

        # Reconstituer les séquences FASTA à partir du fichier Clustal aligné
        merge_sequences(cleaned_clustal_file, merged_fasta_file, transcript_id, imgt_headers)

def remove_header_lines(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        lines = infile.readlines()
        # Skip the first four lines
        for line in lines[3:]:
            outfile.write(line)

def format_fasta_sequence(sequence, line_length=80):
    """
    Formate une séquence en majuscules avec des retours à la ligne tous les `line_length` caractères.
    """
    sequence = sequence.upper()
    return '\n'.join(sequence[i:i + line_length] for i in range(0, len(sequence), line_length))

def align_to_reference(input_file, output_file, reference_id, headers):
    try:
        alignment = AlignIO.read(input_file, "clustal")
        reference_seq = None
        for record in alignment:
            if record.id == reference_id:
                reference_seq = record
                break

        if not reference_seq:
            raise ValueError("Reference sequence not found in the alignment")

        gap_positions = [pos for pos, char in enumerate(str(reference_seq.seq)) if char == "-"]

        with open(output_file, "w") as fasta_file:
            for record in alignment:
                sequence = list(str(record.seq))

                for gap in gap_positions:
                    sequence.insert(gap, "-")

                while sequence[-1] == "-":
                    sequence.pop()

                formatted_sequence = format_fasta_sequence(''.join(sequence))
                header = headers.get(record.id, record.id)
                fasta_file.write(f">{header}\n{formatted_sequence}\n")

    except Exception as e:
        print(f"Erreur lors de la lecture du fichier Clustal: {e}")

def merge_sequences(input_file, output_file, reference_id, headers):
    try:
        alignment = AlignIO.read(input_file, "clustal")
        reference_seq = None
        for record in alignment:
            if record.id == reference_id:
                reference_seq = record
                break

        if not reference_seq:
            raise ValueError("Reference sequence not found in the alignment")

        with open(output_file, "w") as fasta_file:
            for record in alignment:
                if record.id == reference_id:
                    continue
                merged_sequence = []
                ref_seq = str(reference_seq.seq)
                other_seq = str(record.seq)

                for ref_aa, other_aa in zip(ref_seq, other_seq):
                    if ref_aa == "-":
                        merged_sequence.append(other_aa if other_aa != "-" else "")
                    else:
                        merged_sequence.append(ref_aa if ref_aa == other_aa or other_aa == "-" else other_aa)

                final_sequence = ''.join(merged_sequence).replace('-', '')
                formatted_sequence = format_fasta_sequence(final_sequence)
                header = headers.get(record.id, record.id)
                fasta_file.write(f">{header}\n{formatted_sequence}\n")

    except Exception as e:
        print(f"Erreur lors de la lecture du fichier Clustal: {e}")

if __name__ == "__main__":
    parameters = {
        "A": ["ENST00000706902", "forward"],
        "B": ["ENST00000481849", "reverse"],
        "C": ["ENST00000466892", "reverse"],
        "DRB1": ["ENST00000696610", "reverse"],
        "DRB3": ["ENST00000383126", "reverse"],
        "DRB4": ["ENST00000411959", "reverse"],
        "DRB5": ["ENST00000374975", "reverse"],
        "DQA1": ["ENST00000422863", "forward"],
        "DQB1": ["ENST00000487676", "reverse"],
        "DPA1": ["ENST00000692443", "reverse"],
        "DPB1": ["ENST00000418931", "forward"]
    }

    powershell_script = r"D:\mafft-win\mafft-signed.ps1"

    prepare_and_run_mafft(parameters, powershell_script)
