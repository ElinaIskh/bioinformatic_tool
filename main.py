from modules.dna_rna_tools.py import (
    is_nucleic_acid, reverse, transcribe, complement, reverse_complement
    )
from modules.fastq_tools.py import count_gc, count_quality, read_fastq, write_fastq


def run_dna_rna_tools(*args: str) -> "bool | str":
    """
    Takes 1 or more sequences, checks if is it nucleic acid and
    does one of the tools: is_nucleic_acid, reverse, transcribe,
    complement, reverse_complement

    Arguments:
    sequences: list[str]
    tool: str

    Returns bool/str
    Raises exception if:
        Number of argumets is less then 2
        Sequence is not an nucleic acid
        Tool is unknown
    """
    # Checks if there are minimum 2 arguments
    if len(args) < 2:
        raise ValueError(
            "Введите хотя бы одну последовательность и одну процедуру")

    sequences = args[:-1]
    tool = args[-1]

    # Tools dictionary
    tools = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }

    # Checks if the tool exists
    if tool not in tools:
        raise ValueError("Неизвестная процедура")

    #  Checks if the sequence is nucleic acid
    if tool != "is_nucleic_acid":
        for seq in sequences:
            if not is_nucleic_acid(seq):
                raise ValueError("Введена некорректная последовательность")

    #  Resuls of the fuction
    results = [tools[tool](seq) for seq in sequences]
    if len(results) == 1:
        return results[0]
    return results


def filter_fastq(
    input_fastq: str,
    gc_bounds: "tuple | int | float" = (0, 100),
    length_bounds: "tuple | int | float" = (0, 2**32),
    quality_threshold: float = 0,
    output_fastq: str = None
) -> dict:
    """
    Filters reads in fastq file, according to a number of
    guanine-cytosine bounds, length of a read and mean quality
    of a read.

    Arguments:
    input_fastq: str
    gc_bounds: tuple/int/float (default (0, 100))
    length_bounds: tuple/int/float (default (0, 2**32))
    quality_threshold: int/float (default 0)
    output_fastq: str

    Returns str
    """
# Reads the fastq file
    seqs = read_fastq(input_fastq)

# Checks the gc_bounds type and defines the range of GC-bounds
    if isinstance(gc_bounds, (int, float)):
        min_gc = 0
        max_gc = gc_bounds
    elif isinstance(gc_bounds, tuple):
        min_gc = gc_bounds[0]
        max_gc = gc_bounds[1]
    else:
        raise ValueError("Incorrect interval for GC-bounds")

# Checks the length_bounds type and defines the range of sequence's length
    if isinstance(length_bounds, (int, float)):
        min_length = 0
        max_length = length_bounds
    elif isinstance(length_bounds, tuple):
        min_length = length_bounds[0]
        max_length = length_bounds[1]
    else:
        raise ValueError("Incorrect interval for length")

# Filters the reads by input parametrs
    filtered_fastq = {}
    for key, value in seqs.items():
        seq, qual = value
        if (
            (min_gc <= count_gc(seq) <= max_gc)
            and (min_length <= len(seq) <= max_length)
            and count_quality(qual) > quality_threshold
        ):
            filtered_fastq[key] = (seq, qual)

# Writes the fastq file
write_fastq(filtered_fastq)
