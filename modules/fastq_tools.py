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
