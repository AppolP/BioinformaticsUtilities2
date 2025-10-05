def is_nucleic_acid(sequences):
    sequences = set().union(*[s.lower() for s in sequences])
    return sequences.issubset(valid_RNA) or sequences.issubset(valid_DNA)


def transcribe(sequences):
    transcribed_seq = []
    NucAc = ""
    for seq in sequences:
        for nb in seq:
            if nb == "t":
                NucAc += "u"
            elif nb == "T":
                NucAc += "U"
            else:
                NucAc += nb
        transcribed_seq.append(NucAc)
        NucAc = ""
    return transcribed_seq


def reverse(sequences):
    reversed_sequence = []
    for seq in sequences:
        reversed_sequence.append(seq[::-1])
    return reversed_sequence


def complement(sequences):
    sequences = []
    comp_seq = ""
    if sequences.issubset(valid_RNA):
        for seq in sequences:
            sequences = "".join([complement_RNA[nb] for nb in seq])
    elif sequences.issubset(valid_DNA):
        for seq in sequences:
            comp_sequence = "".join([complement_DNA[nb] for nb in seq])
    return comp_sequence


def reverse_complement(sequences):
    return reverse(complement(sequences))


