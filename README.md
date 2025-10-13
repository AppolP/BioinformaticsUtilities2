# Bioinformatics Utilities

A comprehensive Python package for molecular biology and genetic data analysis, providing essential tools for sequence manipulation, reading bioinformatic filed and FASTQ filtering.

## 🧬 DNA/RNA Sequence Tools

### `run_dna_rna_tools()`

Processes DNA and RNA sequences with multiple transformation options. Automatically validates nucleotide composition for all operations.

#### Input

- DNA or RNA sequences with specified transformation type

#### Available Operations

|Operation|Description|
|---|---|
|**`is_nucleic_acid`**|Validates nucleotide composition of sequences|
|**`transcribe`**|Converts DNA sequences to RNA (T → U)|
|**`reverse`**|Reverses sequence direction (5' → 3')|
|**`complement`**|Generates complementary sequences|
|**`reverse_complement`**|Creates reverse complementary sequences|

#### Features

- ✅ Automatic nucleotide validation
    
- ✅ Support for both DNA and RNA sequences
    
- ✅ Comprehensive error handling
    
- ✅ Batch processing capabilities
    

#### Dependencies

- Requires `modules.py` package
    

---

## 🔍 FASTQ Filter

### `filter_fastq()`

Filters sequencing reads based on user-defined quality metrics and criteria.

#### Input Format

```python
{'read_name': ('sequence_string', 'quality_string'),
 'read_1': ('ATCGATCG', 'IIIIIIII'),
 'read_2': ('GCTAGCTA', 'JJJJJJJJ')}
```
#### Filtering Parameters

|Parameter|Type|Default|Description|
|---|---|---|---|
|**`measure_gc`**|`float` (0-100)|`100`|GC content percentage threshold|
|**`measure_quality`**|`int` (Phred score)|`0`|Minimum quality score threshold|
|**`length_bounds`**|`tuple` (min, max)|All reads|Length-based filtering range|

#### Key Metrics

- **GC Content**: Ratio of G and C nucleotides to total sequence length
    
- **Quality Scores**: Phred quality score interpretation
    
- **Length Filtering**: Flexible range-based selection
    

#### Dependencies

- Requires `filter.py` package


---

## 📁 FASTA File Processor

### `convert_multiline_fasta_to_oneline()`

Converts multi-line FASTA files to single-line sequence format. Essential for downstream analysis tools that require single-line sequences.

#### Parameters

|Parameter|Type|Required|Description|
|---|---|---|---|
|**`input_fasta`**|`str`|✅|Path to input FASTA file|
|**`output_fasta`**|`str`|❌ Optional|Path to output FASTA file|

#### Features

- 🧬 Supports DNA, RNA, and protein sequences
    
- 🔄 Converts multi-line sequences to single-line format
    
- ⚡ Optional output file parameter
    
- 🛡️ Automatic validation of FASTA format
    
- 💾 Creates clean, analysis-ready FASTA files
    

#### Example

python

# Basic usage
convert_multiline_fasta_to_oneline("input.fasta", "output.fasta")

# Without output file (creates modified input file)
convert_multiline_fasta_to_oneline("sequences.fasta")

---

## 🔬 BLAST Results Parser

### `parse_blast_output()`

Extracts and processes BLAST analysis results to identify best matches for antibiotic resistance gene analysis in _E. coli_.

#### Parameters

|Parameter|Type|Description|
|---|---|---|
|**`input_file`**|`str`|Path to BLAST output text file|
|**`output_file`**|`str`|Path for saving extracted protein names|

#### Functionality

- 📊 Parses BLAST txt output files
    
- 🔍 Extracts best matches from "Sequences producing significant alignments" sections
    
- 📝 Retrieves first entry from Description column for each QUERY
    
- 📋 Saves alphabetically sorted protein names in single-column format
    

#### Use Case

- **Antibiotic Resistance Analysis**: Identify flanking genes around antibiotic resistance genes in pathogenic _E. coli_
    
- **High-Throughput Processing**: Handles multiple sequence queries in single BLAST run
    
- **Data Organization**: Creates clean, sorted output for further analysis
    

#### Example Input

text

QUERY: gene_flanking_1
Sequences producing significant alignments:
Description          Score    E-value
Protein_A            250      1e-65
Protein_B            200      1e-50

QUERY: gene_flanking_2  
Description          Score    E-value
Protein_C            300      1e-80

---

## 🧫 GBK Gene Extractor

### `select_genes_from_gbk_to_fasta()`

Extracts protein sequences of genes flanking target genes of interest from GenBank files for BLAST analysis.

#### Parameters

|Parameter|Type|Default|Description|
|---|---|---|---|
|**`input_gbk`**|`str`|-|Path to input GBK file|
|**`genes`**|`str` or `list`|-|Target gene(s) of interest|
|**`n_before`**|`int`|`1`|Number of genes to extract before target|
|**`n_after`**|`int`|`1`|Number of genes to extract after target|
|**`output_fasta`**|`str`|-|Output FASTA file path|

#### Features

- 🧬 Extracts protein translations from GBK annotations
    
- 📍 Identifies flanking genes relative to target genes
    
- 🚫 Excludes target genes themselves from output
    
- 🎯 Configurable flanking region size
    
- 🧪 Designed for antibiotic resistance mechanism analysis
    

#### Workflow

1. **Input**: GBK file with _E. coli_ genome annotation
    
2. **Processing**: Locates target genes and extracts n neighboring genes
    
3. **Output**: FASTA file with protein sequences ready for BLAST
    

#### Example Usage

python

# Extract 2 genes before and 1 gene after antibiotic resistance genes
select_genes_from_gbk_to_fasta(
    input_gbk="ecoli_annotation.gbk",
    genes=["blaTEM", "tetA"],
    n_before=2,
    n_after=1,
    output_fasta="flanking_genes.fasta"
)

---

## 🚀 Quick Start

python

from bioinformatics_utils import (
    run_dna_rna_tools, 
    filter_fastq,
    convert_multiline_fasta_to_oneline,
    parse_blast_output,
    select_genes_from_gbk_to_fasta
)

# Sequence transformation
result = run_dna_rna_tools(sequences, operation='reverse_complement')

# FASTQ filtering
filtered_reads = filter_fastq(
    reads_dict, 
    measure_gc=40.0,
    measure_quality=20,
    length_bounds=(50, 150)
)

# FASTA format conversion
convert_multiline_fasta_to_oneline("multiline.fasta", "oneline.fasta")

# BLAST results analysis
parse_blast_output("blast_results.txt", "best_matches.txt")

# Extract flanking genes for resistance analysis
select_genes_from_gbk_to_fasta(
    "ecoli.gbk",
    genes=["resistance_gene"],
    n_before=2,
    n_after=2,
    output_fasta="flanking_sequences.fasta"
)

---

## 💡 Use Cases

- **Sequence Analysis**: Transform and validate DNA/RNA sequences
    
- **Quality Control**: Filter low-quality sequencing reads
    
- **Data Preprocessing**: Prepare FASTA files for analysis
    
- **Antibiotic Resistance Research**: Analyze flanking genes in pathogenic bacteria
    
- **Genome Annotation Processing**: Extract specific gene regions from GBK files
    
- **BLAST Pipeline**: Process and analyze BLAST search results
    
- **Educational Tools**: Teach bioinformatics concepts and workflows
    

---

_Designed for bioinformaticians and molecular biologists working with genetic sequence data, antibiotic resistance research, and genomic analysis._
