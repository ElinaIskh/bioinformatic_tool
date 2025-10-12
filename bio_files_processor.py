import os


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str = None
) -> str:
    """
    Converts the multiline fasta file to oneline fasta file.
    Makes 'output_fasta' file, if there are no such file.

    Arguments:
    input_fasta: str
    output_fasta: str

    Returns str
    """

# Makes 'output_fasta' file, if there are no such file
    if output_fasta is None:
        output_fasta = input_fasta.split(".")[0] + "_output.fasta"

    id = ""
    seq = ""

    with open(input_fasta) as input_file, open(
        output_fasta, "w") as output_file:
        for line in input_file:
            line = line.strip()
            if line.startswith(">"):
                if id != "":
                    output_file.write(id + "\n")
                    output_file.write(seq + "\n")
                id = line
                seq = ""
            else:
                seq += line

# Adds the last pair of id and seq
        if id!= "":
            output_file.write(id + "\n")
            output_file.write(seq + "\n")


def parse_blast_output(
    input_file: str, output_file:str
) -> str:
    """
    Reads txt file (BLAST results), collects all proteins names and sorts
    them in alphabetic order.

    Arguments:
    input_file: str
    output_file: str

    Returns str
    """

    proteins = []

    with open(input_file, "r") as input:
        lines = input.readlines()

    for i in range(len(lines)):
        line = lines[i].strip()
        if line.startswith("Sequences producing significant alignments:"):
            header_id = i + 2  #skips the header
            if header_index < len(lines):
                first_line = lines[header_id].strip()
                protein_name = first_line.split("...")[0].strip()

        if protein_name:
            protein.append(protein_name)

    proteins.sort()

    with open(output_file, "w") as output:
        for protein in proteins:
            output.write(f"{protein}\n")
