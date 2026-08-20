import re
from Bio.Seq import Seq


VALID_BASES = set("ATGC")

COMPLEMENT = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G"
}


def clean_sequence(raw):
    """Uppercase and remove FASTA headers, whitespace, numbers and punctuation."""
    seq = raw.strip().upper()

    lines = seq.splitlines()

    # Remove FASTA header lines
    lines = [line for line in lines if not line.startswith(">")]

    seq = "".join(lines)

    # Keep only letters
    seq = re.sub(r"[^A-Z]", "", seq)

    return seq


def validate_sequence(seq):
    """Return invalid characters in the sequence."""
    invalid = sorted(set(seq) - VALID_BASES)
    return invalid


def complementary_strand(seq):
    """Return the complementary DNA strand."""
    return "".join(COMPLEMENT[base] for base in seq)


def reverse_complement(seq):
    """Return the reverse complement."""
    return complementary_strand(seq)[::-1]


def base_composition(seq):
    """Return A/T/G/C percentages."""
    n = len(seq)

    if not n:
        return {}

    return {
        base: round(seq.count(base) / n * 100, 2)
        for base in "ATGC"
    }


def to_rna(seq):
    """Transcribe DNA to RNA."""
    return seq.replace("T", "U")


def translate_dna(seq):
    """Translate DNA into protein using frame 1."""
    coding_len = len(seq) - (len(seq) % 3)

    protein = str(
        Seq(seq[:coding_len]).translate(to_stop=False)
    )

    return protein


def compare_sequences(seq1, seq2):

#Compares two DNA sequences position by position

    max_length = max(len(seq1), len(seq2))

    padded_seq1 = seq1.ljust(max_length, "-")
    padded_seq2 = seq2.ljust(max_length, "-")

    matches = 0
    mismatches = 0
    gaps = 0

    comparison_line = []

    for base1, base2 in zip(padded_seq1, padded_seq2):

        if base1 == "-" or base2 == "-":

            gaps += 1

            comparison_line.append(" ")

        elif base1 == base2:

            matches += 1

            comparison_line.append("|")

        else:

            mismatches += 1

            comparison_line.append("*")

    aligned_length = max_length

    identity = 0

    if aligned_length > 0:

        identity = round(
            matches / aligned_length * 100,
            2
        )

    return {
        "seq1": seq1,
        "seq2": seq2,
        "aligned_seq1": padded_seq1,
        "aligned_seq2": padded_seq2,
        "comparison_line": "".join(
            comparison_line
        ),
        "length_1": len(seq1),
        "length_2": len(seq2),
        "aligned_length": aligned_length,
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "identity": identity,
        "alignment_score": (
            matches
            - mismatches
            - gaps
        )
    }