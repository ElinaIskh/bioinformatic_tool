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


def filter_fastq(
    seqs: dict,
    gc_bounds: "tuple | int | float" = (0, 100),
    length_bounds: "tuple | int | float" = (0, 2**32),
    quality_threshold: float = 0
) -> dict:
    """
    Filters reads in fastq file, according to a number of
    guanine-cytosine bounds, length of a read and mean
    quality of a read.

    Arguments:
    seqs: dict[str, tuple[str, str]]
    gc_bounds: tuple/int/float (default (0, 100))
    length_bounds: tuple/int/float (default (0, 2**32))
    quality_threshold: int/float (default 0)

    Returns dict[str, tuple[str, str]]
    """


    if isinstance(gc_bounds, (int, float)):
        min_gc = 0
        max_gc = gc_bounds
    elif isinstance(gc_bounds, tuple):
        min_gc = gc_bounds[0]
        max_gc = gc_bounds[1]
    else:
        raise ValueError("Incorrect interval for GC-bounds")

    if isinstance(length_bounds, (int, float)):
        min_length = 0
        max_length = length_bounds
    elif isinstance(length_bounds, tuple):
        min_length = length_bounds[0]
        max_length = length_bounds[1]
    else:
        raise ValueError("Incorrect interval for length")

    filtered_fastq = {}
    for key, value in seqs.items():
        seq, qual = value
        if (
            (min_gc <= count_gc(seq) <= max_gc)
            and (min_length <= len(seq) <= max_length)
            and count_quality(qual) > quality_threshold
        ):
            filtered_fastq[key] = (seq, qual)

    return filtered_fastq
