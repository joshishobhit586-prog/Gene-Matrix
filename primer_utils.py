COMPLEMENT = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G"
}


def reverse_complement(seq):
    return "".join(COMPLEMENT[base] for base in reversed(seq))


def gc_percent(seq):
    if not seq:
        return 0.0

    gc = seq.count("G") + seq.count("C")

    return round(gc / len(seq) * 100, 2)


def melting_temp(seq):
    """Approximate primer melting temperature (degrees C).

    Uses the simple Wallace rule for very short primers and the
    GC%-based formula commonly used for standard-length PCR primers
    (~18-25 nt). This is an estimate for screening candidates, not a
    substitute for a full nearest-neighbor thermodynamic calculation.
    """
    n = len(seq)

    if n == 0:
        return 0.0

    a = seq.count("A")
    t = seq.count("T")
    g = seq.count("G")
    c = seq.count("C")

    if n < 14:
        return float(2 * (a + t) + 4 * (g + c))

    return round(64.9 + 41 * (g + c - 16.4) / n, 2)


def has_run(seq, length=4):
    """True if the sequence contains a run of the same base
    `length` or longer (a common cause of mispriming)."""
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
    """A 3' G or C ('GC clamp') improves primer binding stability."""
    return seq[-1] in ("G", "C")


def score_primer(seq, target_gc=(40, 60), target_tm=(55, 65)):
    """Lower is better. Penalizes GC% and Tm outside the target
    ranges, missing GC clamp, and homopolymer runs."""
    gc = gc_percent(seq)
    tm = melting_temp(seq)

    penalty = 0.0

    if gc < target_gc[0]:
        penalty += (target_gc[0] - gc)
    elif gc > target_gc[1]:
        penalty += (gc - target_gc[1])

    if tm < target_tm[0]:
        penalty += (target_tm[0] - tm)
    elif tm > target_tm[1]:
        penalty += (tm - target_tm[1])

    if not gc_clamp_ok(seq):
        penalty += 2

    if has_run(seq, 4):
        penalty += 5

    return penalty


def find_candidate_primers(seq, start, end, length_range=(18, 25),
                            target_gc=(40, 60), target_tm=(55, 65)):
    """Slide a window of varying length across seq[start:end] and
    return a list of candidate primer dicts (sorted best-first)."""
    candidates = []

    min_len, max_len = length_range

    for primer_len in range(min_len, max_len + 1):
        for pos in range(start, end - primer_len + 1):
            primer = seq[pos:pos + primer_len]

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
                "penalty": score_primer(primer, target_gc, target_tm)
            })

    candidates.sort(key=lambda c: c["penalty"])

    return candidates


def design_primer_pairs(seq, length_range=(18, 25),
                         product_range=(150, 500),
                         target_gc=(40, 60), target_tm=(55, 65),
                         search_window=200, max_pairs=5):
    """Design forward/reverse PCR primer pair candidates for `seq`.

    Forward primers are searched near the 5' end and reverse primers
    near the 3' end, so the resulting amplicon spans (most of) the
    input sequence. Returns up to `max_pairs` pairs, best first.
    """
    n = len(seq)

    fwd_end = min(n, search_window)
    rev_start = max(0, n - search_window)

    fwd_candidates = find_candidate_primers(
        seq, 0, fwd_end, length_range, target_gc, target_tm
    )

    rev_region_candidates = find_candidate_primers(
        seq, rev_start, n, length_range, target_gc, target_tm
    )

    pairs = []

    # Cap the search space so this stays fast for longer sequences
    for fwd in fwd_candidates[:60]:
        for rev_region in rev_region_candidates[:60]:

            product_size = rev_region["end"] - fwd["start"]

            if product_size < product_range[0] or product_size > product_range[1]:
                continue

            if rev_region["start"] < fwd["end"]:
                continue

            reverse_primer_seq = reverse_complement(rev_region["seq"])

            tm_diff = abs(fwd["tm"] - rev_region["tm"])

            pair_penalty = fwd["penalty"] + rev_region["penalty"] + tm_diff

            pairs.append({
                "forward_seq": fwd["seq"],
                "forward_start": fwd["start"] + 1,  # 1-based for display
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

    pairs.sort(key=lambda p: p["penalty"])

    return pairs[:max_pairs]
