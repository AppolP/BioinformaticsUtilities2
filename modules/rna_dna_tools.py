from typing import Union

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
