# Bioinformatic tool

A small collection of Python utilities for basic manipulation of DNA/RNA sequences and simple FASTQ read filtering.  

Author: Elina Iskhakova, *Institute of Cytology, Saint-Petersburg, Russia*

---

## Project structure

```
bioinformatic_tool/
│
├── README.md
├── main.py                      # Main script (entry point)
└── modules/
    ├── dna_rna_tools.py         # DNA/RNA sequence tools
    └── fastq_tools.py           # FASTQ utilities (GC content, quality)

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

## FASTQ filtering  `filter_fastq`

Filters reads by GC content, read length, and mean base quality.
```
filter_fastq(
    seqs,
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold=0
)
```

**Arguments:**

- `seqs`: `dict[str, tuple[str, str]]` — dictionary of reads, e.g.
    ```
    {
    "read1": ("ATGCGT", "IIIIII"),
    "read2": ("ATTTGC", "II!!II")
}

- `gc_bounds`: allowed GC content range (%)
- `length_bounds`: allowed sequence length range
- `quality_threshold`: minimum mean quality

**Returns:** dictionary of filtered reads.

**Example:**

```
>>> reads = {
...     "r1": ("ATGC", "IIII"),
...     "r2": ("ATAT", "!!!!")
... }
>>> filter_fastq(reads, gc_bounds=(40, 60), quality_threshold=30)
{'r1': ('ATGC', 'IIII')}
```

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

- `count_gc(seq)` — calculates GC percentage
- `count_quality(qual)` — calculates mean quality score

