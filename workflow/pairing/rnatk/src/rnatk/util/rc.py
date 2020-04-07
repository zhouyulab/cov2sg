def reverse_complement(dna):
    """Reverse complement of given DNA"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join([complement[base] for base in reversed(dna)])


def reverse(dna):
    return "".join(reversed(dna))
