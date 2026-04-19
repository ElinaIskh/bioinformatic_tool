import pytest
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from seq_tool import DNASequence, RNASequence, parse_bounds, filter_by_gc, write_fastq


class TestSequences:
    def test_dna_validity(self):
        dna = DNASequence("ATGC")
        assert dna.is_valid() is True

    def test_transcription(self):
        dna = DNASequence("ATGC")
        assert str(dna.transcribe()) == "AUGC"

    def test_parse_bounds_error(self):
        with pytest.raises(ValueError):
            parse_bounds([1, 2, 3])


class TestFiltering:
    def test_gc_filter(self):
        rec = SeqRecord(Seq("GGGG"), id="test")  # GC 100%
        assert filter_by_gc(rec, 0, 50) is False
        assert filter_by_gc(rec, 80, 100) is True

    def test_sequence_indexing(self):
        dna = DNASequence("ATGC")
        assert str(dna[0:2]) == "AT"


class TestIO:
    def test_file_creation(self, tmp_path):
        os.chdir(tmp_path)
        rec = SeqRecord(Seq("ATGC"), id="test", letter_annotations={
                        "phred_quality": [40, 40, 40, 40]})
        path = write_fastq([rec], "test.fastq")
        assert os.path.exists(path)

    def test_log_exists(self):
        assert os.path.exists("tool_log.log")

    def test_reverse_complement(self):
        dna = DNASequence("ATGC")
        assert str(dna.reverse_complement()) == "GCAT"
