from modules.rna_dna_tools import is_nucleic_acid, transcribe, reverse, complement, reverse_complement
from modules.filter import is_in_bounds, is_qualified
complement_RNA: dict
complement_RNA = {
    "a": "u",
    "A": "U",
    "u": "a",
    "U": "A",
    "g": "c",
    "G": "C",
    "c": "g",
    "C": "G",
}
complement_DNA: dict
complement_DNA = {
    "a": "t",
    "A": "T",
    "t": "a",
    "T": "A",
    "g": "c",
    "G": "C",
    "c": "g",
    "C": "G",
}
valid_RNA: dict
valid_RNA = {"a", "u", "g", "c", "A", "U", "G", "C"}
valid_DNA: dict
valid_DNA = {"a", "t", "g", "c", "A", "T", "G", "C"}

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
    Raises the error in case of incorrect sequences
    """
    sequences = args[:-1]
    action = str(args[-1])
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
    else:
        answer = available_functions[action](sequences)
        if len(answer) == 1:
            return str(answer[0])
        else:
            return answer
        

def filter_fastq(
    seqs: dict[str: str], 
    gc_bounds: Union[tuple[int, int], int, float] = (0, 100), 
    length_bounds: Union[tuple[int, int], int, float] = (0, 2**32), 
    quality_threshold: int = 0
) -> dict[str, str]:
    """
    Filters fastq sequences by gives parameters
    Arguments:
    seqs: dict 
    gc_bounds: tuples / int / float
    length_bounds: tuples / int / float
    quality_threshold: int
    Returns dict
    Raises the error in case reads are not correct nucleic acids
    """
    for key in seqs:
        if not is_nucleic_acid(key):
            return 'Error: reads are not nucleic acids'
    filtered_seqs = {}
    if isinstance(length_bounds, (int, float)):
        len_left_bound, len_right_bound = (0, length_bounds)
    else:
        len_left_bound, len_right_bound = length_bounds
    for seq in seqs:
        if is_in_bounds(seqs[seq][0], gc_bounds, 'G', 'C') and is_qualified(seqs[seq][1], quality_threshold) and len_left_bound <= len(seqs[seq][0]) <= len_right_bound:
            filtered_seqs[seq] = seqs[seq]
    return filtered_seqs
          
    