
def is_nucleic_acid(sequences: Union[str, list]) -> bool:
    sequences = set().union(*[s.lower() for s in sequences])
    return sequences.issubset(valid_rna) or sequences.issubset(valid_dna)


def transcribe(sequences: Union[str, list]) -> list:
    return [''.join([transcribed_dna[nucleotide] for nucleotide]) for seq in sequences]


def reverse(sequences: Union[str, list]) -> list:
    reversed_sequences = [seq[::-1] for seq in sequences]
    return reversed_sequences


def complement(sequences: Union[str, list]) -> list:
    if sequences.issubset(valid_rna):
        complement_sequences = [''.join([complement_rna[nb] for nb in seq]) for seq in sequences]
    elif sequences.issubset(valid_dna):
        complement_sequences = [''.join([complement_dna[nb] for nb in seq]) for seq in sequences]
    return complement_sequences


def reverse_complement(sequences: Union[str, list]) -> list: 
    return reverse(complement(sequences))


