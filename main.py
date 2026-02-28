import os
from abc import ABC, abstractmethod
from typing import Union, List
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction


class BiologicalSequence(ABC):
    """
    Abstract base class for biological sequences.
    
    Provides a common interface for sequence length, indexing, 
    string representation, and alphabet validation.
    """

    def __init__(self, seq: str):
        self.seq = seq

    def __len__(self) -> int:
        return len(self.seq)

    def __getitem__(self, index: Union[int, slice]) -> "BiologicalSequence":
        return self.__class__(self.seq[index])

    def __str__(self) -> str:
        return self.seq

    @abstractmethod
    def is_valid(self) -> bool:
        """
        Checks if the sequence consists of valid characters for its type.
        """
        pass


class NucleicAcidSequence(BiologicalSequence):
    """
    Represents nucleic acid sequences (DNA and RNA).
    
    Implements core operations like complementation and reversal
    based on a specific complement map.
    """

    def __init__(self, seq: str):
        super().__init__(seq)

    def is_valid(self) -> bool:
        """Verify if all characters belong to the defined nucleic acid alphabet."""
        return set(self.seq.upper()) <= self.alphabet

    def complement(self) -> "NucleicAcidSequence":
        """
        Return the complement sequence using the class's complement map.
        Raises NotImplementedError if the map is not defined.
        """
        if not hasattr(self, "complement_map"):
            raise NotImplementedError(
                "Method .compliment() is not known for this class")
    
        new_seq = self.seq.translate(self.complement_map)
        return self.__class__(new_seq)

    def reverse(self) -> "NucleicAcidSequence":
        """Return the reversed sequence"""
        return self.__class__(self.seq[::-1])

    def reverse_complement(self) -> "NucleicAcidSequence":
        """Return the reverse-complement of the sequence."""
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    """Class for DNA sequences."""

    alphabet = {"A", "T", "G", "C", "a", "t", "g", "c"}
    complement_map = str.maketrans("ATGCatgc", "TACGtacg")

    def __init__(self, seq):
        super().__init__(seq)

    def transcribe(self) -> "RNASequence":
        """Transcribe the DNA sequence into an RNA sequence."""
        new_seq = self.seq.replace("T", "U").replace("t", "u")
        return RNASequence(new_seq)


class RNASequence(NucleicAcidSequence):
    """Class for RNA sequences."""

    alphabet = {"A", "U", "G", "C", "a", "u", "g", "c"}
    complement_map = str.maketrans("AUGCaugc", "UACGuacg")

    def __init__(self, seq):
        super().__init__(seq)


class AminoAcidSequence(BiologicalSequence):
    """Class for Amino acid sequences."""

    alphabet = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy")

    def __init__(self, seq):
        super().__init__(seq)

    def is_valid(self) -> bool:
        """Verify if all characters belong to the standard amino acid alphabet."""
        return set(self.seq) <= self.alphabet

    def count_aminoacid(self, aminoacid: str) -> int:
        """Count the occurrences of a specific amino acid in the sequence."""
        return self.seq.count(aminoacid)


def filter_fastq(
    input_fastq: str,
    gc_bounds: Union[tuple, int, float] = (0, 100),
    length_bounds: Union[tuple, int, float] = (0, 2**32),
    quality_threshold: float = 0,
) -> List[SeqRecord]:
    """
    Filters reads from a FASTQ file, according to a number of
    guanine-cytosine bounds, length of a read and mean quality
    of a read. 

    Arguments:
    input_fastq: str
    gc_bounds: tuple/int/float (default (0, 100))
    length_bounds: tuple/int/float (default (0, 2**32))
    quality_threshold: int/float (default 0)
    
    Returns a list of SeqRecord objects.
    """
    min_gc, max_gc = parse_bounds(gc_bounds, "GC")
    min_len, max_len = parse_bounds(length_bounds, "length")

    filtered_fastq = []
    for record in SeqIO.parse(input_fastq, "fastq"):
        if (
            filter_by_gc(record, min_gc, max_gc)
            and filter_by_length(record, min_len, max_len)
            and filter_by_quality(record, quality_threshold)
        ):
            filtered_fastq.append(record)
    
    return filtered_fastq


def parse_bounds(bounds, name = "parameter"):
    """
    Parses bounds (tuple/int/float) and returns (min, max). Paramenter: GC or length.
    Raises ValueError if bounds are invalid.
    """
    if isinstance(bounds, (int, float)):
        return 0, bounds
    elif isinstance(bounds, tuple) and len(bounds) == 2:
        return bounds
    else:
        raise ValueError(f"Incorrect interval for {name}-bounds: {bounds}")


def filter_by_gc(record, min_gc: float, max_gc: float) -> bool:
    """Checks if GC-content of the sequence is within the bounds."""
    gc = gc_fraction(record.seq) * 100 
    return min_gc <= gc <= max_gc


def filter_by_length(record, min_len: int, max_len: int) -> bool:
    """Checks if sequence length is within the bounds."""
    return min_len <= len(record.seq) <= max_len


def filter_by_quality(record, threshold: float) -> bool:
    """Checks if mean quality exceeds the threshold."""
    qualities = record.letter_annotations["phred_quality"]
    mean_qual = sum(qualities) / len(qualities)
    return mean_qual >= threshold


def write_fastq(records: List[SeqRecord], output_filename: str) -> None:
    """
    Writes a list of SeqRecord objects to a FASTQ file in the 'filtered' directory.
    """
    if not os.path.exists("filtered"):
        os.makedirs("filtered")
    output_path = os.path.join("filtered", output_filename)
    
    count = SeqIO.write(records, output_path, "fastq")
    print(f"Successfully wrote {count} records to {output_path}")

