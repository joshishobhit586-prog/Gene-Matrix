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
    (~18-25 nt).
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
    """True if the sequence contains a run of the same base."""

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
    """A 3' G or C improves primer binding stability."""
    return bool(seq) and seq[-1] in ("G", "C")


def score_primer(seq, target_gc=(40, 60), target_tm=(55, 65)):
    """Lower score is better."""

    gc = gc_percent(seq)
    tm = melting_temp(seq)
    penalty = 0.0

    if gc < target_gc[0]:
        penalty += target_gc[0] - gc
    elif gc > target_gc[1]:
        penalty += gc - target_gc[1]

    if tm < target_tm[0]:
        penalty += target_tm[0] - tm
    elif tm > target_tm[1]:
        penalty += tm - target_tm[1]

    if not gc_clamp_ok(seq):
        penalty += 2

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
    """Design forward/reverse PCR primer pair candidates."""

    n = len(seq)

    fwd_end = min(n, search_window)
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

    # Cap the search space so this stays fast
    for fwd in fwd_candidates[:60]:

        for rev_region in rev_region_candidates[:60]:

            product_size = (
                rev_region["end"] - fwd["start"]
            )

            if (
                product_size < product_range[0]
                or product_size > product_range[1]
            ):
                continue

            if rev_region["start"] < fwd["end"]:
                continue

            reverse_primer_seq = reverse_complement(
                rev_region["seq"]
            )

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

    pairs.sort(key=lambda p: p["penalty"])

    selected = []
    used_forward = set()
    used_reverse = set()

    for p in pairs:
        f = p["forward_seq"]
        r = p["reverse_seq"]

    # Never reuse the same forward primer
        if f in used_forward:
            continue

    # Never reuse the same reverse primer
        if r in used_reverse:
            continue

        selected.append(p)
        used_forward.add(f)
        used_reverse.add(r)

        if len(selected) >= max_pairs:
            break

    return selected