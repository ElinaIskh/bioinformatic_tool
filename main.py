from .modules.dna_rna_tools import (
    is_nucleic_acid, reverse, transcribe, complement, reverse_complement
    )
from .modules.fastq_tools import count_gc, count_quality, read_fastq, write_fastq


def run_dna_rna_tools(*args: str) -> "bool | str":
    """
    Takes 1 or more sequences, checks if is it nucleic acid and
    does one of the tools: is_nucleic_acid, reverse, transcribe,
    complement, reverse_complement.
    This functions checks if there are at least 2 arguments, 
    if the tool exists and if the sequence is nucleic acid.

    Arguments:
    sequences: list[str]
    tool: str

    Returns bool/str
    Raises exception if:
        Number of argumets is less then 2
        Sequence is not an nucleic acid
        Tool is unknown
    """

    if len(args) < 2:
        raise ValueError(
            "Write at least one sequence and one tool")

    *sequences, tool = args

    # Tools dictionary
    tools = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }

    if tool not in tools:
        raise ValueError("Unknown tool")

    if tool != "is_nucleic_acid":
        for seq in sequences:
            if not is_nucleic_acid(seq):
                raise ValueError("Incorrect sequence")

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

    seqs = read_fastq(input_fastq)

    min_gc, max_gc = parse_bounds(gc_bounds, "GC")
    min_len, max_len = parse_bounds(length_bounds, "length")

    filtered_fastq = {}
    for key, (seq, qual) in seqs.items():
        if (
            filter_by_gc(seq, min_gc, max_gc)
            and filter_by_length(seq, min_len, max_len)
            and filter_by_quality(qual, quality_threshold)
        ):
            filtered_fastq[key] = (seq, qual)
    return filtered_fastq


# Writes the fastq file
filtered_fastq = filter_fastq('bioinformatic_tool\example_data\example_fastq.fastq', quality_threshold=15, gc_bounds=50, length_bounds=500)

write_fastq(filtered_fastq, "my_output.fastq")
