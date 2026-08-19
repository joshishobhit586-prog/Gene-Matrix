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


def global_alignment(seq1, seq2):
    """
    Perform Needleman-Wunsch global alignment.

    Scoring:
    Match = +1
    Mismatch = -1
    Gap = -1
    """

    match_score = 1
    mismatch_score = -1
    gap_score = -1

    rows = len(seq1) + 1
    cols = len(seq2) + 1

    score_matrix = [
        [0 for _ in range(cols)]
        for _ in range(rows)
    ]

    # Initialize first column
    for i in range(1, rows):
        score_matrix[i][0] = i * gap_score

    # Initialize first row
    for j in range(1, cols):
        score_matrix[0][j] = j * gap_score

    # Fill matrix
    for i in range(1, rows):

        for j in range(1, cols):

            if seq1[i - 1] == seq2[j - 1]:
                diagonal = (
                    score_matrix[i - 1][j - 1]
                    + match_score
                )
            else:
                diagonal = (
                    score_matrix[i - 1][j - 1]
                    + mismatch_score
                )

            up = (
                score_matrix[i - 1][j]
                + gap_score
            )

            left = (
                score_matrix[i][j - 1]
                + gap_score
            )

            score_matrix[i][j] = max(
                diagonal,
                up,
                left
            )

    # Traceback
    aligned_seq1 = []
    aligned_seq2 = []

    i = len(seq1)
    j = len(seq2)

    while i > 0 or j > 0:

        if i > 0 and j > 0:

            if seq1[i - 1] == seq2[j - 1]:
                diagonal_score = (
                    score_matrix[i - 1][j - 1]
                    + match_score
                )
            else:
                diagonal_score = (
                    score_matrix[i - 1][j - 1]
                    + mismatch_score
                )

            if score_matrix[i][j] == diagonal_score:

                aligned_seq1.append(
                    seq1[i - 1]
                )

                aligned_seq2.append(
                    seq2[j - 1]
                )

                i -= 1
                j -= 1

                continue

        if (
            i > 0
            and score_matrix[i][j]
            == score_matrix[i - 1][j] + gap_score
        ):

            aligned_seq1.append(
                seq1[i - 1]
            )

            aligned_seq2.append("-")

            i -= 1

        else:

            aligned_seq1.append("-")

            aligned_seq2.append(
                seq2[j - 1]
            )

            j -= 1

    aligned_seq1.reverse()
    aligned_seq2.reverse()

    return (
        "".join(aligned_seq1),
        "".join(aligned_seq2)
    )


def compare_sequences(seq1, seq2):
    """
    Compare two DNA sequences using global alignment.
    """

    aligned_seq1, aligned_seq2 = global_alignment(
        seq1,
        seq2
    )

    matches = 0
    mismatches = 0
    gaps = 0

    comparison_line = []

    for base1, base2 in zip(
        aligned_seq1,
        aligned_seq2
    ):

        if base1 == base2 and base1 != "-":

            matches += 1

            comparison_line.append("|")

        elif base1 == "-" or base2 == "-":

            gaps += 1

            comparison_line.append(" ")

        else:

            mismatches += 1

            comparison_line.append("*")

    aligned_length = len(aligned_seq1)

    identity = 0

    if aligned_length > 0:

        identity = round(
            matches / aligned_length * 100,
            2
        )

    return {
        "seq1": seq1,
        "seq2": seq2,
        "aligned_seq1": aligned_seq1,
        "aligned_seq2": aligned_seq2,
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