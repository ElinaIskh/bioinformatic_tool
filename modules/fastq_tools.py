import os

def count_gc(seq: str) -> float:
    """
    Counts number of guanine-cytosine bounds in sequence.

    Arguments:
    seq: str

    Returns float
    """
    g_number = seq.count("G")
    c_number = seq.count("C")
    nucleotide_number = len(seq)
    return round((g_number + c_number)/nucleotide_number * 100, 2)


def count_quality(qual: str) -> float:
    """
    Counts mean quality of sequence.

    Arguments:
    seq: str

    Returns float
    """
    char_quality = 0
    for char in qual:
        char_quality += ord(char) - 33
    return round(char_quality/len(qual), 2)


def read_fastq(input_fastq: str) -> dict:
    """
    Reads the input fastq file and turns it into the dictionary:
    key: ID of sequence; value: (sequence, quality).

    Arguments:
    input_fastq: str

    Returns dict[str, tuple[str, str]]
    """
    seqs = {}
    with open(input_fastq) as fastq_file:
        while True:
            id = fastq_file.readline()
            if not id:
                break
            seq = fastq_file.readline()
            plus = fastq_file.readline()
            qual = fastq_file.readline()
            seqs[id.strip()] = (seq.strip(), qual.strip())
    return seqs


def write_fastq(filtered_fastq: dict, output_fastq: str) -> str:
    """
    Writes the filtered dictionary of sequences into the fastq file.
    Makes 'filtered' directory, if there are no such directory.
    Does not rewrite output_fastq file, if it already exists.

    Arguments:
    filtered_fastq: dict
    output_fastq: str

    Returns str
    Raises exception if output_fastq file already exists.
    """
    if not os.path.exists("filtered"):
        os.mkdir("filtered")
    output_path = os.path.join("filtered", output_fastq)
    if os.path.exists(output_path):
        print("File output_fastq already exists")
        return

    with open(output_path, "w") as fastq_file:
        for id, (seq, qual) in filtered_fastq.items():
            fastq_file.write(f"{id}\n{seq}\n + \n{qual}\n")
