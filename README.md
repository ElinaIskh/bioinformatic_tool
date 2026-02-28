# Bioinformatic tool

A collection of Python utilities for 

- manipulation of DNA/RNA sequences 
- FASTQ read filtering.
- FASTA filtering
- BLAST readingbioinformatic_tool/

Author: Elina Iskhakova, *Institute of Cytology, Saint-Petersburg, Russia*

---

## Project structure

```
bioinformatic_tool/
│
├── README.md
├── main.py                      # OOP Classes & FASTQ Filter (Core Logic)
└── bio_files_processor.py       # Fasta and BLAST file utilities

```

## Object-Oriented Sequence Manipulation

The tool provides a robust class hierarchy for handling biological data as objects rather than simple strings.

**Supported Classes:**

- `DNASequence`: Supports `complement()`, `reverse()`, `reverse_complement()`, and `transcribe()`.
- `RNASequence`: Supports `complement()`, `reverse()`, and `reverse_complement()`
- `AminoAcidSequence`: Supports `count_aminoacid(symbol)`.

**Key Capabilities:**

- Validation: Automatically check if sequences match their biological alphabet.
- Slicing: Slice objects (e.g., dna[1:5]) and receive a new object of the same class.
- Polymorphism: Native handling of complementation rules without manual type checking.

Example:
``` python
from main import DNASequence

dna = DNASequence("ATGC")
print(dna.reverse_complement())  # Returns a DNASequence object: GCAT
rna = dna.transcribe()           # Returns an RNASequence object: AUGC
```

## FASTQ Filtering (Biopython Integrated)

Efficiently filter sequencing reads using the `Bio.SeqIO` engine. Filters are applied to GC content, read length, and mean Phred quality.

Function: `filter_fastq`


``` python
filter_fastq(
    input_fastq: str,
    gc_bounds: Union[tuple, int, float] = (0, 100),
    length_bounds: Union[tuple, int, float] = (0, 2**32),
    quality_threshold: float = 0
) -> List[SeqRecord]
```

**Arguments:**

- `input_fastq`: path to the input fastq file
- `gc_bounds`: allowed GC content range (%)
- `length_bounds`: allowed sequence length range
- `quality_threshold`: minimum mean quality

**Returns:** fastq file

**Example:**

```python 
from main import filter_fastq, write_fastq

records = filter_fastq("data.fastq", quality_threshold=30)
write_fastq(records, "filtered_data.fastq")  # Saves to 'filtered/' directory
```

## Bio file processing — modules/bio_files_processor.py
### Converting multiline FASTA to oneline

``` python
convert_multiline_fasta_to_oneline("input.fasta", "output.fasta")
```

- input_fasta — input file

- output_fasta — output file (optional; default adds _output.fasta)

Converts sequences in FASTA to one-line per sequence format.

### Parsing BLAST output

``` python
parse_blast_output("blast_results.txt", "proteins_sorted.txt")
```

Reads a BLAST output text file, extracts protein names, sorts them alphabetically, and writes them to output_file.

---

## Installation & Requirements
This tool requires Python 3.8+ and the Biopython library.
```bash
pip install biopython
```

Import the required classes or functions directly from main.py or bio_files_processor.py to integrate them into your pipeline.
```python
from main import DNASequence, filter_fastq
```
