from modules.rna_dna_tools import is_nucleic_acid, transcribe, reverse, complement, reverse_complement
from modules.filter import is_in_bounds, is_qualified
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
valid_RNA = {"a", "u", "g", "c", "A", "U", "G", "C"}
valid_DNA = {"a", "t", "g", "c", "A", "T", "G", "C"}

def run_dna_rna_tools(*args):
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
        

def filter_fastq(seqs, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0):
    filtered_seqs = {}
    if isinstance(length_bounds, (int, float)):
        len_left_bound, len_right_bound = (0, length_bounds)
    else:
        len_left_bound, len_right_bound = length_bounds
    for seq in seqs:
        if is_in_bounds(seqs[seq][0], gc_bounds, 'G', 'C') and is_qualified(seqs[seq][1], quality_threshold) and len_left_bound <= len(seqs[seq][0]) <= len_right_bound:
            filtered_seqs[seq] = seqs[seq]
    return filtered_seqs
          
    