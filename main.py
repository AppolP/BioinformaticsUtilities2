from modules.rna_dna_tools import is_nucleic_acid, transcribe, reverse, complement, reverse_complement
from modules.filter import is_in_bounds, is_qualified
from typing import Union
import os

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
    
    *sequences, action= args
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
    output_fastq: str = 'filtered/filtered_fastq.fastq'
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
    
    if not os.path.exists('filtered'):
        os.makedirs('filtered')
    if os.path.exists('filtered/output_fastq'):
        return 'File already exists!'
    
    if isinstance(length_bounds, (int, float)):
        len_left_bound, len_right_bound = (0, length_bounds)
    else:
        len_left_bound, len_right_bound = length_bounds    
        
    with open(input_fastq, 'r') as raw_fastq, open(output_fastq, 'w'):
        output_fastq.write(f"{input_fastq = }{'\n'}")
        
        skip_comma = True
        while raw_fastq.readline().strip() != '':
            name_line = raw_fastq.readline().strip(':')[0]
            sequence = raw_fastq.readline().strip()
            comment = raw_fastq.readline()
            quality = raw_fastq.readline().strip()
            
            if not is_nucleic_acid(sequence):
                return 'Error: reads are not nucleic acids'
            name = name_line[1:] if name_line.startswith('@') else name_line
            if (is_in_bounds(sequence, gc_bounds, 'G', 'C') and 
                is_qualified(quality, quality_threshold) and 
                len_left_bound <= len(sequence) <= len_right_bound):
                if not skip_comma:
                    output_fastq.write(',\n')
                output_fastq.write(f"'{name}': ('{sequence}', '{quality}')")
                skip_comma = False
                
        output_fastq.write("/n}}")
       
    # for key in seqs:
    #     if not is_nucleic_acid(key):
    #         return 'Error: reads are not nucleic acids'
    # filtered_seqs = {}
    
    # for seq in seqs:
    #     if is_in_bounds(seqs[seq][0], gc_bounds, 'G', 'C') and is_qualified(seqs[seq][1], quality_threshold) and len_left_bound <= len(seqs[seq][0]) <= len_right_bound:
    #         filtered_seqs[seq] = seqs[seq]
    # return filtered_seqs
          
    