from abc import ABC, abstractmethod
from typing import Union


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




