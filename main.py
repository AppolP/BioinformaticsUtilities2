from modules.filter import is_in_bounds, is_qualified
from typing import Union
import os


valid_rna: dict
valid_rna = {"a", "u", "g", "c", "A", "U", "G", "C"}

valid_dna: dict
valid_dna = {"a", "t", "g", "c", "A", "T", "G", "C"}

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

complement_rna: dict
complement_rna = {
    "a": "u",
    "A": "U",
    "u": "a",
    "U": "A",
    "g": "c",
    "G": "C",
    "c": "g",
    "C": "G",
}

complement_dna: dict
complement_dna = {
    "a": "t",
    "A": "T",
    "t": "a",
    "T": "A",
    "g": "c",
    "G": "C",
    "c": "g",
    "C": "G",
}


def is_nucleic_acid(sequences: Union[str, list]) -> bool:
    """
    Сhecks whether the sequence contains only nucleic acids

    Arguments:
    sequence: str, list, sequence to process

    Returns True or False
    """

    sequences = set().union(*[s.lower() for s in sequences])
    return sequences.issubset(valid_rna) or sequences.issubset(valid_dna)


def transcribe(sequences: Union[str, list]) -> list:
    """
    Exchanges T (thymin) with U (uracil) the corresponding amino acid in the RNA

    Arguments:
    sequence: str, list, sequence or sequences to process

    Returns list with transcribed sequence or sequences
    """

    return [
        "".join([transcribed_dna[nucleotide] for nucleotide in seq])
        for seq in sequences
    ]


def reverse(sequences: Union[str, list]) -> list:
    """
    Reads the sequence backwards

    Arguments:
    sequence: str, list, sequence or sequences to process

    Returns list with reversed sequence or sequences
    """

    reversed_sequences = [seq[::-1] for seq in sequences]
    return reversed_sequences


def complement(sequences: Union[str, list]) -> list:
    """
    Exchanges each nucleic acid with its complemented pair

    Arguments:
    sequence: str, list, sequence or sequences to process

    Returns list with complemented sequence or sequences
    """

    if sequences.issubset(valid_rna):
        complement_sequences = [
            "".join([complement_rna[nb] for nb in seq]) for seq in sequences
        ]
    elif sequences.issubset(valid_dna):
        complement_sequences = [
            "".join([complement_dna[nb] for nb in seq]) for seq in sequences
        ]
    return complement_sequences


def reverse_complement(sequences: Union[str, list]) -> list:
    """
    Exchanges each nucleic acid with its complemented pair and reads it backwards

    Arguments:
    sequence: str, list, sequence or sequences to process

    Returns list with complemented and turned in the opposite direction sequence or sequences
    """

    return reverse(complement(sequences))

def is_in_bounds(sequence: str, range: Union[tuple[int, int], int, float], *args: str) -> bool:
    '''
    Check if sequence length is in given bounds
    
    Arguments:
    sequence: str, sequence to process
    range: tuple, int or float, a range the sequence belongs to or a value the sequence is below 
    args: str, substrings to be counted and matched to the specified range
    '''
    if isinstance(range, (int, float)):
        low_bound, upper_bound = 0, range
    else:
        low_bound, upper_bound = range
    args_count = 0
    for arg in args:
        args_count += sequence.count(arg)  
    args_percentage = args_count/len(sequence)*100
    return low_bound <= args_percentage <= upper_bound

def is_qualified(seq_quality: str, quality_threshold: int) -> bool:
    '''
    Check if sequence quality under the given threshold
    
    Arguments: 
    seq_quality: str, quality sequence of the read sequence 
    quality_threshold: int, the value bordering acceptable quality value
    '''
    q_score = sum([ord(el)-33 for el in seq_quality])  
    return 10**(-q_score/10) >= quality_threshold

def run_dna_rna_tools(*args: str) -> Union[str, list[str]]:
    """
    Processes rna or dna according to given order by last string in agrs

    Arguments:
    args[:-1]: str, sequences to process
    agrs[-1]: str, what to do with a sequence
    Returns str in case of a single sequence or list in case of 2 and more sequences
    Possuble args[-1]:
    is_nucleic_acid - check for correct nucleic acid
    transcribe - transform DNA to RNA
    reverse - write the sequence in the opposite direction
    complement - write complement sequence
    reverse_complement - combine reverse and complement transformations

    Raises an error in case of incorrect sequences
    """

    *sequences, action = args
    available_functions = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }
    flag = is_nucleic_acid(sequences)
    if not flag:
        return flag
    answer = available_functions[action](sequences)
    if len(answer) == 1:
        return str(answer[0])
    return answer


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
    input_fastq: str
    gc_bounds: tuples / int / float
    length_bounds: tuples / int / float
    quality_threshold: int
    output_fastq: str

    Returns filtered fastq sequences in a user-specified directory. Saves results in the current directory by default
    Raises the error in case reads are not correct nucleic acids
    """

    if not os.path.exists("filtered"):
        os.makedirs("filtered")
    if os.path.exists("filtered/output_fastq"):
        return "File already exists!"

    if isinstance(length_bounds, (int, float)):
        len_left_bound, len_right_bound = (0, length_bounds)
    else:
        len_left_bound, len_right_bound = length_bounds

    with open(input_fastq, "r") as raw_fastq, open(output_fastq, "w") as output_fastq:
        while True:
            name_line = raw_fastq.readline().strip()
            if not name_line:
                break

            sequence = raw_fastq.readline().strip()
            plus_line = raw_fastq.readline().strip()
            quality = raw_fastq.readline().strip()

            if not name_line.startswith("@"):
                continue
            if not is_nucleic_acid(sequence):
                print("Error: reads are not nucleic acids")
                return None

            name = name_line[1:]

            if (
                is_in_bounds(sequence, gc_bounds, "G", "C")
                and (is_qualified(quality, quality_threshold))
                and (len_left_bound <= len(sequence) <= len_right_bound)
            ):

                output_fastq.write(f"{name_line}\n")
                output_fastq.write(f"{sequence}\n")
                output_fastq.write(f"{plus_line}\n")
                output_fastq.write(f"{quality}\n")
