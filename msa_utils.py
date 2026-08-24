import re


VALID_BASES = set("ATGC")


def parse_multi_fasta(text):
    """Parse pasted text into a list of (name, sequence) tuples.

    Supports standard multi-FASTA (">name" headers) as well as plain
    text where each non-empty line is treated as its own sequence.
    """
    text = text.strip()

    if not text:
        return []

    if ">" in text:
        entries = []
        blocks = text.split(">")

        for block in blocks:
            block = block.strip()

            if not block:
                continue

            lines = block.splitlines()
            name = lines[0].strip() or f"Seq{len(entries) + 1}"
            seq = "".join(lines[1:])
            seq = re.sub(r"[^A-Za-z]", "", seq).upper()

            if seq:
                entries.append((name, seq))

        return entries

    # No FASTA headers: one sequence per non-empty line
    entries = []

    for i, line in enumerate(text.splitlines(), start=1):
        seq = re.sub(r"[^A-Za-z]", "", line).upper()

        if seq:
            entries.append((f"Seq{i}", seq))

    return entries


def validate_sequences(entries):
    """Return a dict of {name: [invalid chars]} for any sequence with
    characters outside A/T/G/C."""
    invalid = {}

    for name, seq in entries:
        bad = sorted(set(seq) - VALID_BASES)

        if bad:
            invalid[name] = bad

    return invalid


def needleman_wunsch(seq1, seq2, match=2, mismatch=-1, gap=-2):
    """Classic global pairwise alignment (Needleman-Wunsch).

    Returns (aligned_seq1, aligned_seq2, score).
    """
    n, m = len(seq1), len(seq2)

    score = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        score[i][0] = i * gap

    for j in range(1, m + 1):
        score[0][j] = j * gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score[i - 1][j - 1] + (
                match if seq1[i - 1] == seq2[j - 1] else mismatch
            )
            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap

            score[i][j] = max(diag, up, left)

    # Traceback
    aligned1 = []
    aligned2 = []

    i, j = n, m

    while i > 0 or j > 0:

        if i > 0 and j > 0:
            current = score[i][j]
            diag = score[i - 1][j - 1] + (
                match if seq1[i - 1] == seq2[j - 1] else mismatch
            )

            if current == diag:
                aligned1.append(seq1[i - 1])
                aligned2.append(seq2[j - 1])
                i -= 1
                j -= 1
                continue

        if i > 0 and score[i][j] == score[i - 1][j] + gap:
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1
            continue

        aligned1.append("-")
        aligned2.append(seq2[j - 1])
        j -= 1

    aligned1.reverse()
    aligned2.reverse()

    return "".join(aligned1), "".join(aligned2), score[n][m]


def _gap_run_slots(aligned_center, n):
    """Given one pairwise alignment row for the center sequence
    (length n before gaps were inserted), return a list of length
    n + 1 with the number of gap characters that fall in each slot
    (slot k = gaps immediately before original center position k;
    slot n = trailing gaps after the last character)."""
    slots = [0] * (n + 1)

    slot = 0
    run = 0

    for ch in aligned_center:
        if ch == "-":
            run += 1
        else:
            slots[slot] = run
            run = 0
            slot += 1

    slots[slot] = run

    return slots


def star_alignment(entries):
    """Simplified progressive ("center-star") multiple sequence
    alignment. Not a full ClustalW implementation, but produces a
    reasonable multiple alignment using only pure Python.

    entries: list of (name, seq)
    Returns (names, aligned_seqs) where aligned_seqs are all the
    same length and share consistent gap columns.
    """
    if len(entries) < 2:
        raise ValueError("At least two sequences are required.")

    names = [name for name, _ in entries]
    seqs = [seq for _, seq in entries]

    # Use the longest sequence as the center - a common, simple
    # heuristic that tends to minimize gaps introduced into it.
    center_idx = max(range(len(seqs)), key=lambda i: len(seqs[i]))
    center_seq = seqs[center_idx]
    n = len(center_seq)

    pairwise = []  # (aligned_center_i, aligned_other_i) for every other sequence

    for i, seq in enumerate(seqs):
        if i == center_idx:
            continue

        aligned_c, aligned_s, _ = needleman_wunsch(center_seq, seq)
        pairwise.append((i, aligned_c, aligned_s))

    # Determine, for every gap "slot" around the center sequence, the
    # maximum number of gap characters used by any pairwise alignment.
    max_gaps = [0] * (n + 1)
    slot_data = {}

    for i, aligned_c, aligned_s in pairwise:
        slots = _gap_run_slots(aligned_c, n)
        slot_data[i] = slots

        for k in range(n + 1):
            if slots[k] > max_gaps[k]:
                max_gaps[k] = slots[k]

    # Build the final center row
    final_center_parts = []

    for k in range(n):
        final_center_parts.append("-" * max_gaps[k])
        final_center_parts.append(center_seq[k])

    final_center_parts.append("-" * max_gaps[n])
    final_center = "".join(final_center_parts)

    # Build the final row for every other sequence
    results = {center_idx: final_center}

    for i, aligned_c, aligned_s in pairwise:
        final_row_parts = []
        run_chars = []
        slot = 0

        for cc, sc in zip(aligned_c, aligned_s):
            if cc == "-":
                run_chars.append(sc)
            else:
                pad = max_gaps[slot] - len(run_chars)
                final_row_parts.append("".join(run_chars) + "-" * pad)
                final_row_parts.append(sc)
                run_chars = []
                slot += 1

        pad = max_gaps[slot] - len(run_chars)
        final_row_parts.append("".join(run_chars) + "-" * pad)

        results[i] = "".join(final_row_parts)

    aligned_seqs = [results[i] for i in range(len(seqs))]

    return names, aligned_seqs


def consensus_sequence(aligned_seqs):
    """Majority-vote consensus across aligned columns. Gaps only
    appear in the consensus if every sequence has a gap in that
    column."""
    if not aligned_seqs:
        return ""

    length = len(aligned_seqs[0])
    consensus = []

    for col in range(length):
        counts = {}

        for seq in aligned_seqs:
            base = seq[col]
            counts[base] = counts.get(base, 0) + 1

        # Prefer a real base over a gap when there's a tie
        best = max(
            counts.items(),
            key=lambda item: (item[1], item[0] != "-")
        )[0]

        consensus.append(best)

    return "".join(consensus)


def conservation_line(aligned_seqs):
    """Per-column conservation symbols, loosely inspired by Clustal's
    output: '*' = fully conserved (all sequences share the same base,
    no gaps), ' ' = otherwise."""
    if not aligned_seqs:
        return ""

    length = len(aligned_seqs[0])
    line = []

    for col in range(length):
        column = {seq[col] for seq in aligned_seqs}

        if len(column) == 1 and "-" not in column:
            line.append("*")
        else:
            line.append(" ")

    return "".join(line)


def percent_identity(seq1, seq2):
    """Percent identity between two equal-length aligned sequences,
    over columns where at least one sequence has a real base."""
    compared = 0
    matches = 0

    for a, b in zip(seq1, seq2):
        if a == "-" and b == "-":
            continue

        compared += 1

        if a == b:
            matches += 1

    if compared == 0:
        return 0.0

    return round(matches / compared * 100, 2)


def pairwise_identity_matrix(names, aligned_seqs):
    """Return an NxN list-of-lists of percent identity values."""
    n = len(names)
    matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(n):
            if i == j:
                matrix[i][j] = 100.0
            else:
                matrix[i][j] = percent_identity(aligned_seqs[i], aligned_seqs[j])

    return matrix
