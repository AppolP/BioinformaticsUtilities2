import os

from typing import Union
from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


class BiologicalSequence(ABC):


    def __init__(self, sequence: str):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)
    
    def __getitem__(self, index: Union[int, slice]):
            if isinstance(index, slice): 
                start, stop, step = index.start, index.stop, index.step
                insequence = BiologicalSequence(self.sequence[start:stop:step])
                return insequence
            else :
                insequence = BiologicalSequence(self.sequence[index])
                return insequence
        
    def __str__(self):
        return self.sequence
    



class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, sequence: str):
        super().__init__(sequence)

    complement_library = {}
    valid_rna: set
    valid_rna = {"a", "u", "g", "c", "A", "U", "G", "C"}

    valid_dna: set
    valid_dna = {"a", "t", "g", "c", "A", "T", "G", "C"}

    def complement(self):
        if not self.complement_library:
            raise NotImplementedError
        return NucleicAcidSequence("".join([self.complement_library[nb] for nb in self.sequence]))

    def reverse(self):
        return NucleicAcidSequence(self.sequence[::-1])

    def reverse_complement(self):
        return NucleicAcidSequence(self.complement()[::-1])
    
    def is_correct(self):
        sequence_set = list(self.sequence)
        return sequence_set.issubset(self.valid_rna) or sequence_set.issubset(self.valid_dna)
    

class DNASequence(NucleicAcidSequence):
    complement_library: dict
    complement_library = {
        "a": "t",
        "A": "T",
        "t": "a",
        "T": "A",
        "g": "c",
        "G": "C",
        "c": "g",
        "C": "G",
    }
    transcribed_dna: dict
    transcribed_dna = {
        "a": "a",
        "A": "A",
        "t": "u",
        "T": "U",
        "g": "g",
        "G": "G",
        "c": "c",
        "C": "C",
    }

    def __init__(self, sequence: str):
        super().__init__(sequence)
    
    def transcribe(self):
        return RNASequence("".join([self.transcribed_dna[nucleotide] for nucleotide in self.sequence]))

class RNASequence(NucleicAcidSequence):
    complement_library: dict
    complement_library = {
        "a": "u",
        "A": "U",
        "u": "a",
        "U": "A",
        "g": "c",
        "G": "C",
        "c": "g",
        "C": "G",
    }

    def __init__(self, sequence: str):
        super().__init__(sequence)
    


class AminoAcidSequence(BiologicalSequence):
    amino_table = {
        'UUU': 'F', 'UUC': 'F',  
        'UUA': 'L', 'UUG': 'L', 
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',  
        'UAU': 'Y', 'UAC': 'Y',  
        'UAA': '*', 'UAG': '*',  
        'UGU': 'C', 'UGC': 'C',  
        'UGA': '*',              
        'UGG': 'W',              
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',  
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',  
        'CAU': 'H', 'CAC': 'H', 
        'CAA': 'Q', 'CAG': 'Q',  
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I',  
        'AUG': 'M',            
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',  
        'AAU': 'N', 'AAC': 'N',  
        'AAA': 'K', 'AAG': 'K',  
        'AGU': 'S', 'AGC': 'S',  
        'AGA': 'R', 'AGG': 'R',  
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',  
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',  
        'GAU': 'D', 'GAC': 'D', 
        'GAA': 'E', 'GAG': 'E',  
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
    
    def is_correct(self):
        if not self.sequence:
            return False 
        valid_amino = set(self.amino_table.values())
        return all(ch.upper() in valid_amino for ch in self.sequence)
    
    def __init__(self, sequence: str):
        super().__init__(sequence)

    def protein_synthesis(self):
        if not self.is_correct():
            raise ValueError("Invalid amino acid sequence")
        
        rna = RNASequence(self.sequence.upper())
        if not rna.is_correct():
            raise ValueError("Invalid RNA sequence")
        
        amino_sequence = []

        for i in range(0, len(rna)-2, 3):
            codon = rna[i:i+3]
            if codon == '*':
                break
            amino_sequence.append(self.amino_table[codon])
        return AminoAcidSequence(''.join(amino_sequence))


def filter_fastq(
    input_fastq: str,
    gc_bounds: Union[tuple[int, int], int, float] = (0, 100),
    length_bounds: Union[tuple[int, int], int, float] = (0, 2**32),
    quality_threshold: int = 0,
    output_fastq: str = "filtered/filtered_fastq.fastq",
) -> None:
    """
    Read and filter fastq sequences in a user-specified dierectory by given parameters

    Arguments:
    input_fastq: str / list 
    gc_bounds: tuples / int / float
    length_bounds: tuples / int / float
    quality_threshold: int
    output_fastq: str

    Returns filtered fastq sequences in a user-specified directory. Saves results in the current directory by default
    Raises the error in case reads are not correct nucleic acids
    """

    if not os.path.exists(output_fastq):
        os.makedirs(os.path.dirname(output_fastq))

    if isinstance(length_bounds, (int, float)):
        len_left_bound, len_right_bound = (0, length_bounds)
    else:
        len_left_bound, len_right_bound = length_bounds

    if isinstance(gc_bounds, (int, float)):
        low_bound, upper_bound = (0, gc_bounds)
    else:
        low_bound, upper_bound = gc_bounds
    
    output_file = output_fastq.split('/')[-1]
    output_pw = output_fastq.split('/')[:-1]

    with open(input_fastq, "r") as raw_fastq, open("/".join(output_pw + [output_file]), "w") as output_fastq:

        sequences = SeqIO.parse(raw_fastq, "fastq")
        filtered = []

        for sequence in sequences:
            seq = str(sequence.seq).upper()
            if not len_left_bound <= len(sequence) <= len_right_bound:
                continue

            gc_content = gc_fraction(seq) * 100
            if not (low_bound <= gc_content <= upper_bound):
                continue
            
            qualities = sequence.letter_annotations["phred_quality"]
            if qualities:
                mean_quality = sum(qualities) / len(qualities)
                if mean_quality < quality_threshold:
                    continue
            elif quality_threshold > 0:
                continue

            filtered.append(sequence)
        
        SeqIO.write(filtered, output_fastq, "fastq")



