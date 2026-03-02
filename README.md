# Bioinformatics Utilities

A comprehensive Python package for molecular biology and genetic data analysis, providing essential tools for sequence manipulation, reading bioinformatic files and FASTQ filtering.  
**All biological sequences are now implemented as classes** (`DNASequence`, `RNASequence`, `AminoAcidSequence`) with built-in methods for reverse complement, transcription, translation, and validity checks.  
The FASTQ filter has been rewritten to use **Biopython** for robust and efficient processing.

---

## Sequence Classes

The package provides object‑oriented wrappers for DNA, RNA, and amino acid sequences.

| Class | Description | Key Methods |
|-------|-------------|-------------|
| `DNASequence(sequence: str)` | DNA sequence (letters A,T,G,C) | `complement()`, `reverse()`, `reverse_complement()`, `transcribe()` |
| `RNASequence(sequence: str)` | RNA sequence (letters A,U,G,C) | `complement()`, `reverse()`, `reverse_complement()` |
| `AminoAcidSequence(sequence: str)` | Protein sequence (one‑letter code) | `protein_synthesis()` (translates RNA into protein), `is_correct()` |

All classes inherit from `BiologicalSequence` and support:
- `len()` – length of the sequence
- indexing and slicing (`seq[2]`, `seq[1:5]`)
- `str(seq)` – returns the raw sequence string
- `is_correct()` – checks if all characters belong to the allowed alphabet

### Example usage

```python

dna = DNASequence("ATGGCC")
print(dna.reverse_complement())          # GGCCAT
print(dna.transcribe())                   # AUGGCC

rna = RNASequence("AUGGCC")
print(rna.complement())                    # UACCGG
print(rna.reverse())                        # CCGGUA

protein = AminoAcidSequence("MVLSPADKT")
print(protein.is_correct())                 # True
```

---

## FASTQ Filter

### `filter_fastq()`

Filters sequencing reads from a FASTQ file using GC content, length, and average quality thresholds. **Now powered by Biopython** for reliable parsing and writing.

#### Input

- Path to a FASTQ file (plain or gzipped).

#### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| **`gc_bounds`** | `tuple(min, max)` or `int/float` | `(0, 100)` | Allowed GC content range (percentage). If a single number is given, it is treated as the upper bound (0–value). |
| **`length_bounds`** | `tuple(min, max)` or `int/float` | `(0, 2**32)` | Allowed sequence length range. If a single number is given, it becomes the upper bound (0–value). |
| **`quality_threshold`** | `int` | `0` | Minimum average Phred quality score required. |
| **`output_fastq`** | `str` | `"filtered/filtered_fastq.fastq"` | Path where filtered reads will be saved. The directory is created automatically if missing. |

#### Returns

- `None` – writes filtered reads directly to the output file.

#### Example

```python

filter_fastq(
    input_fastq="raw_reads.fastq",
    gc_bounds=(30, 70),
    length_bounds=(100, 500),
    quality_threshold=20,
    output_fastq="filtered/high_quality.fastq"
)
```

---

## FASTA File Processor

### `convert_multiline_fasta_to_oneline()`

Converts multi-line FASTA files to single-line sequence format. Essential for downstream analysis tools that require single-line sequences.

#### Parameters

|Parameter|Type|Required|Description|
|---|---|---|---|
|**`input_fasta`**|`str`|✅|Path to input FASTA file|
|**`output_fasta`**|`str`|❌ Optional|Path to output FASTA file|

#### Features

- Supports DNA, RNA, and protein sequences
- Converts multi-line sequences to single-line format
- Optional output file parameter (if omitted, the input file is overwritten)
- Automatic validation of FASTA format
- Creates clean, analysis-ready FASTA files

#### Example

```python
from bio_files_processor import convert_multiline_fasta_to_oneline

# Without output file (modifies input file in place)
convert_multiline_fasta_to_oneline("multiline.fasta")

# With explicit output
convert_multiline_fasta_to_oneline("multiline.fasta", "oneline.fasta")
```
bio_files_processor
---

## BLAST Results Parser

### `parse_blast_output()`

Extracts and processes BLAST analysis results to identify best matches for antibiotic resistance gene analysis in *E. coli*.

#### Parameters

|Parameter|Type|Description|
|---|---|---|
|**`input_file`**|`str`|Path to BLAST output text file|
|**`output_file`**|`str`|Path for saving extracted protein names|

#### Functionality

- Parses BLAST txt output files
- Extracts best matches from "Sequences producing significant alignments" sections
- Retrieves first entry from Description column for each QUERY
- Saves alphabetically sorted protein names in single-column format

#### Example

```python
from bio_files_processor import parse_blast_output

parse_blast_output("blast_results.txt", "best_matches.txt")
```

---

## Quick Start

```python

# --- Sequence manipulation ---
dna = DNASequence("ATGGCC")
print(dna.reverse_complement())          # GGCCAT
print(dna.transcribe())                   # AUGGCC

rna = RNASequence("AUGGCC")
print(rna.complement())                    # UACCGG

# --- FASTQ filtering ---
filter_fastq(
    input_fastq="raw.fastq",
    gc_bounds=50,
    length_bounds=(100, 300),
    quality_threshold=25,
    output_fastq="filtered/pass.fastq"
)

# --- FASTA conversion ---
convert_multiline_fasta_to_oneline("multiline.fasta", "oneline.fasta")

# --- BLAST parsing ---
parse_blast_output("blast.txt", "best_hits.txt")

```

---

## Dependencies

- **Python** ≥ 3.7
- **Biopython** (for FASTQ filtering)

All other modules are part of the Python standard library.

To install the required package:

```bash
pip install biopython
```

Or use the provided `requirements.txt`:

```
biopython
```

---

## Use Cases

- **Sequence Analysis**: Transform and validate DNA/RNA sequences using an intuitive class‑based interface.
- **Quality Control**: Filter low‑quality sequencing reads with flexible GC, length, and quality criteria.
- **Data Preprocessing**: Prepare FASTA files for analysis (multi‑line → single‑line).
- **Antibiotic Resistance Research**: Analyze flanking genes in pathogenic bacteria using GenBank files.
- **BLAST Pipeline**: Process and extract best matches from BLAST search results.
- **Educational Tools**: Teach bioinformatics concepts and workflows.

---

*Designed for bioinformaticians and molecular biologists working with genetic sequence data, antibiotic resistance research, and genomic analysis.*