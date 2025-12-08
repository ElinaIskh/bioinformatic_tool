# Bioinformatic tool

A collection of Python utilities for 

- manipulation of DNA/RNA sequences 
- FASTQ read filtering.
- FASTA filtering
- BLAST protein sorting

Author: Elina Iskhakova, *Institute of Cytology, Saint-Petersburg, Russia*

---

## Project structure

```
bioinformatic_tool/
│
├── README.md
├── main.py                      # Main script (entry point)
├── bio_files_processor.py       # Fasta and BLAST file utilities
└── modules/
    ├── dna_rna_tools.py         # DNA/RNA sequence tools
    └── fastq_tools.py           # FASTQ processing and quality control

```

## DNA/RNA utilities  `run_dna_rna_tools`

This function applies various operations to nucleic acid sequences.
`run_dna_rna_tools(seq1, seq2, ..., tool)`

**Arguments:**

- `seq1, seq2, ...` — one or more nucleotide sequences
- `tool` — name of the operation to perform  
    (`is_nucleic_acid`, `reverse`, `transcribe`, `complement`, `reverse_complement`)

**Example:**

```
>>> run_dna_rna_tools("AUGC", "reverse") 
"CGUA"  

>>> run_dna_rna_tools("ATGC", "reverse_complement")
"GCAT"  

>>> run_dna_rna_tools("ATUGC", "is_nucleic_acid") 
False
```

## Reading and writing FASTQ

```
from modules.fastq_tools import read_fastq, write_fastq

seqs = read_fastq("input.fastq")
write_fastq(seqs, "output.fastq")  # writes into 'filtered/' directory
```

Notes:

- read_fastq(input_fastq: str) -> dict[str, tuple[str, str]] — reads a FASTQ file into a dictionary with keys as sequence IDs and values as (sequence, quality) tuples.

- write_fastq(filtered_fastq: dict, output_fastq: str) -> str — writes sequences to a FASTQ file in a filtered/ directory. Does not overwrite existing files.


### FASTQ filtering  `filter_fastq`

Filters reads by GC content, read length, and mean base quality.
```
filter_fastq(
    input_fastq,
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold=0,
    output_fastq
)
```

**Arguments:**

- `input_fastq`: name of input fastq file
- `gc_bounds`: allowed GC content range (%)
- `length_bounds`: allowed sequence length range
- `quality_threshold`: minimum mean quality
- `output_fastq`: name of output fastq file

**Returns:** fastq file

**Example:**

```
>>> filter_fastq('example_fastq.fastq', gc_bounds=(40, 60), quality_threshold=30, 'out_fastq.fastq')
```

## Bio file processing — modules/bio_files_processor.py
### Converting multiline FASTA to oneline

```
convert_multiline_fasta_to_oneline("input.fasta", "output.fasta")
```

- input_fasta — input file

- output_fasta — output file (optional; default adds _output.fasta)

Converts sequences in FASTA to one-line per sequence format.

### Parsing BLAST output

```
parse_blast_output("blast_results.txt", "proteins_sorted.txt")
```

Reads a BLAST output text file, extracts protein names, sorts them alphabetically, and writes them to output_file.

## Modules overview

### `modules/dna_rna_tools.py`

Contains:

- `is_nucleic_acid(seq)` — checks if the input is a valid nucleic acid
- `reverse(seq)` — returns the reversed sequence
- `transcribe(seq)` — transcribes DNA to RNA 
- `complement(seq)` — returns the complementary strand
- `reverse_complement(seq)` — returns the reverse complement strand

### `modules/fastq_tools.py`

Contains:

- `read_fastq(input_fastq)` — read FASTQ file into a dictionary
- `write_fastq(filtered_fastq, output_fastq)` — write FASTQ dictionary to file
- `count_gc(seq)` — calculates GC percentage
- `count_quality(qual)` — calculates mean quality score
- Filtering helpers: `filter_by_gc()`, `filter_by_length()`, `filter_by_quality()`
- `parse_bounds()` — validate and parse numeric bounds



