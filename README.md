# DNA/RNA Sequence Parser

The Python script in this repo parses FASTA files to determine sequence types (DNA/RNA), identify open reading frames (ORFs), and compute sequence statistics.

For each sequence, the script identifies its type (DNA/RNA).
It finds Open Reading Frames (ORFs) based on start codons (ATG/AUG) and stop codons (TAA, TAG, TGA for DNA, and UAA, UAG, UGA for RNA).
Computes statistics like sequence length, base counts, and number of ORFs. Before we dive into the scripts, let's get some understanding of the domain.

## What is Bioinformatics?

Bioinformatics is an interdisciplinary field of science that develops methods and software tools for understanding biological data, especially when the data sets are large and complex.

### FASTA Format

FASTA format is a text-based format for representing either nucleotide sequences or amino acid (protein) sequences. Nucleotides or amino acids are represented using single-letter codes. Link to full article [here](https://en.wikipedia.org/wiki/FASTA_format)

## Format/Syntax:
A sequence begins with a greater-than character (">") followed by a name and/or description of the sequence (all in a single line).
The lines immediately following the description line are the sequence representation, with one letter per amino acid or nucleic acid, and are typically no more than 80 characters in length. (80 character length is not considered in the Python script)

## How to Run the Python Script

`python dna_parser.py`

## How to Run the Unittest Script

`python test_dna_parser.py`

## Assumptions

DNA Sequences: Contains A, C, G, T, and ambiguous bases (N, R, Y, etc.).

RNA Sequences: Contains A, C, G, U, and ambiguous bases.

NOTE: referred to the wikipedia page on [FASTA](https://en.wikipedia.org/wiki/FASTA_format) for determining the ambiguous characters. Assumed ambiguous characters are: N, R, Y, K, M, S, W, B, D, H, V, X, Z

## Repository Structure
```bash
/
│
├── test_files/                  # Folder containing test files
│   ├── test_part1_in.fasta      # Example test file 1
│   ├── test_part2_in.fasta      # Example test file 2
│   ├── test_part3_in.fasta      # Example test file 3
│
├── dna_parser.py                # Main Python script for parsing DNA/RNA sequences
├── test_dna_parser.py           # Unit tests for the dna_parser.py script
└── README.md                    # Documentation for how to run and use the script
```

