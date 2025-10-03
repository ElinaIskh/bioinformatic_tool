def is_nucleic_acid(seq: str) -> bool:
    """
    Checks if the imput is nucleic acid or not

    Arguments:
    seq: str

    Returns bool
    Raises exception if seq is not an nucleic acid
    """
    seq_upper = seq.upper()
    if "T" in seq_upper and "U" in seq_upper:
        return False
    dna = set("ATGC")
    rna = set("AUGC")
    return (set(seq_upper) <= dna) or (set(seq_upper) <= rna)


def reverse(seq: str) -> str:
    """
    Makes a reverse sequence of DNA or RNA

    Arguments:
    seq: str

    Returns str
    """
    return seq[::-1]


def transcribe(seq: str) -> str:
    """
    Makes a transcribed sequence of DNA or RNA

    Arguments:
    seq: str

    Returns str
    """
    seq = seq.replace("t", "u")
    seq = seq.replace("T", "U")
    return seq


def complement(seq: str) -> str:
    """
    Makes a complement sequence of DNA or RNA

    Arguments:
    seq: str

    Returns str
    """
    table = str.maketrans("ATGCatgc", "TACGtacg")
    return seq.translate(table)


def reverse_complement(seq: str) -> str:
    """
    Makes a reverse complement sequence of DNA or RNA

    Arguments:
    seq: str

    Returns str
    """
    return reverse(complement(seq))
