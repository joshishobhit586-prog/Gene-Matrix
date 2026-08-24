COMPLEMENT = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G"
}


def complement(seq):
    """Return the direct complementary DNA sequence."""
    return "".join(COMPLEMENT[base] for base in seq)


def reverse_complement(seq):
    """Return the reverse complementary DNA sequence."""
    return "".join(COMPLEMENT[base] for base in reversed(seq))


def gc_percent(seq):
    """Calculate GC percentage."""
    if not seq:
        return 0.0

    gc = seq.count("G") + seq.count("C")
    return round(gc / len(seq) * 100, 2)


def melting_temp(seq):
    """Approximate primer melting temperature in degrees Celsius."""

    n = len(seq)

    if n == 0:
        return 0.0

    a = seq.count("A")
    t = seq.count("T")
    g = seq.count("G")
    c = seq.count("C")

    # Wallace rule for short primers
    if n < 14:
        return float(2 * (a + t) + 4 * (g + c))

    # Approximate formula for standard-length primers
    return round(64.9 + 41 * (g + c - 16.4) / n, 2)


def has_run(seq, length=4):
    """Check for a run of the same nucleotide."""

    if not seq:
        return False

    run = 1

    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            run += 1

            if run >= length:
                return True
        else:
            run = 1

    return False


def gc_clamp_ok(seq):
    """Check whether the primer ends with G or C."""

    return bool(seq) and seq[-1] in ("G", "C")


def score_primer(seq, target_gc=(40, 60), target_tm=(55, 65)):
    """Score a primer. Lower score is better."""

    gc = gc_percent(seq)
    tm = melting_temp(seq)

    penalty = 0.0

    # GC percentage penalty
    if gc < target_gc[0]:
        penalty += target_gc[0] - gc

    elif gc > target_gc[1]:
        penalty += gc - target_gc[1]

    # Melting temperature penalty
    if tm < target_tm[0]:
        penalty += target_tm[0] - tm

    elif tm > target_tm[1]:
        penalty += tm - target_tm[1]

    # Prefer a GC clamp
    if not gc_clamp_ok(seq):
        penalty += 2

    # Penalize homopolymer runs
    if has_run(seq, 4):
        penalty += 5

    return penalty


def find_candidate_primers(
    seq,
    start,
    end,
    length_range=(18, 25),
    target_gc=(40, 60),
    target_tm=(55, 65)
):
    """Find and rank candidate primers."""

    candidates = []

    min_len, max_len = length_range

    for primer_len in range(min_len, max_len + 1):

        for pos in range(start, end - primer_len + 1):

            # Extract region from original input sequence
            template_region = seq[pos:pos + primer_len]

            # Primer is the direct complement of the input sequence region
            # Example: GGTGC -> CCACG
            primer = complement(template_region)

            if has_run(primer, 5):
                continue

            candidates.append({
                "seq": primer,
                "start": pos,
                "end": pos + primer_len,
                "length": primer_len,
                "gc": gc_percent(primer),
                "tm": melting_temp(primer),
                "gc_clamp": gc_clamp_ok(primer),
                "penalty": score_primer(
                    primer,
                    target_gc,
                    target_tm
                )
            })

    candidates.sort(key=lambda c: c["penalty"])

    return candidates


def design_primer_pairs(
    seq,
    length_range=(18, 25),
    product_range=(150, 500),
    target_gc=(40, 60),
    target_tm=(55, 65),
    search_window=200,
    max_pairs=5,
    max_reuse_per_primer=1
):
    """Design and rank forward/reverse primer pairs."""

    n = len(seq)

    # Need enough sequence to create a primer
    if n < length_range[0]:
        return []

    # Search forward primers near the beginning
    fwd_end = min(n, search_window)

    # Search reverse primers near the end
    rev_start = max(0, n - search_window)

    fwd_candidates = find_candidate_primers(
        seq,
        0,
        fwd_end,
        length_range,
        target_gc,
        target_tm
    )

    rev_region_candidates = find_candidate_primers(
        seq,
        rev_start,
        n,
        length_range,
        target_gc,
        target_tm
    )

    pairs = []

    # Cap search space for performance
    for fwd in fwd_candidates[:60]:

        for rev_region in rev_region_candidates[:60]:

            product_size = (
                rev_region["end"] - fwd["start"]
            )

            # Product must be within the requested size range
            if (
                product_size < product_range[0]
                or product_size > product_range[1]
            ):
                continue

            # Prevent overlapping forward and reverse regions
            if rev_region["start"] < fwd["end"]:
                continue

            # The reverse candidate is already the complement
            # of its corresponding input sequence region
            reverse_primer_seq = rev_region["seq"]

            # Compare melting temperatures
            tm_diff = abs(
                fwd["tm"] - rev_region["tm"]
            )

            pair_penalty = (
                fwd["penalty"]
                + rev_region["penalty"]
                + tm_diff
            )

            pairs.append({
                "forward_seq": fwd["seq"],
                "forward_start": fwd["start"] + 1,
                "forward_end": fwd["end"],
                "forward_gc": fwd["gc"],
                "forward_tm": fwd["tm"],

                "reverse_seq": reverse_primer_seq,
                "reverse_start": rev_region["start"] + 1,
                "reverse_end": rev_region["end"],
                "reverse_gc": rev_region["gc"],
                "reverse_tm": rev_region["tm"],

                "product_size": product_size,
                "tm_diff": round(tm_diff, 2),
                "penalty": round(pair_penalty, 2)
            })

    # Best-scoring pairs first
    pairs.sort(key=lambda p: p["penalty"])

    # Select unique primer pairs
    selected = []

    used_forward = {}
    used_reverse = {}

    for p in pairs:

        f = p["forward_seq"]
        r = p["reverse_seq"]

        # Prevent excessive reuse of the same forward primer
        if used_forward.get(f, 0) >= max_reuse_per_primer:
            continue

        # Prevent excessive reuse of the same reverse primer
        if used_reverse.get(r, 0) >= max_reuse_per_primer:
            continue

        selected.append(p)

        used_forward[f] = used_forward.get(f, 0) + 1
        used_reverse[r] = used_reverse.get(r, 0) + 1

        if len(selected) >= max_pairs:
            break

    return selected