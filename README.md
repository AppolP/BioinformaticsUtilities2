Bioinformatics Utilities
This package provides two essential functions for molecular biology and genetic data analysis.

1. DNA/RNA Sequence Tools (run_dna_rna_tools)
This function processes DNA and RNA sequences with multiple transformation options. It automatically validates nucleotide composition for all operations.

Input: DNA or RNA sequences with specified transformation type:

is_nucleic_acid - validates nucleotide composition

transcribe - converts DNA to RNA sequences

reverse - reverses sequence direction

complement - generates complementary sequences

reverse_complement - creates reverse complementary sequences

Dependencies: Requires the modules.py package.

2. FASTQ Filter (filter_fastq)
This function filters sequencing reads based on user-defined criteria.

Input: Dictionary in format {'name': ('read', 'quality')} with optional parameters:

measure_gc - GC content percentage (ratio of G and C nucleotides to total sequence length). Default: 100 (sequences with 100% GC content)

measure_quality - Phred quality score threshold. Default: 0

length_bounds - length-based read filtering. Default: all reads pass

Dependencies: Requires the filter.py package.
