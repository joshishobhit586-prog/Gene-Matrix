import streamlit as st

from dna_utils import clean_sequence, validate_sequence
from primer_utils import design_primer_pairs


st.title("PCR Primer Designer")

st.write(
    "Enter a DNA sequence to design forward and reverse PCR primer "
    "candidates for amplifying it."
)

st.caption(
    "Primers are ranked by GC content, melting temperature (Tm), a 3' "
    "GC clamp, and matched Tm between the pair. Tm is estimated with a "
    "standard GC%-based approximation, not a full nearest-neighbor "
    "thermodynamic model — always confirm critical primers with a "
    "dedicated tool before ordering."
)


example_seq = (
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCGTGGGGAAGCTGCATGT"
    "GACCGATTAGCGCTAGCTAGCTAGCATCGATCGTAGCTAGCTAGCATCGATCGATCGTAG"
    "CATCGATCGTAGCTAGCATCGATCGTAGCTAGCTAGCATCGATCGATCGTAGCTAGCATC"
    "GTAGCTAGCATCGATCGATCGATCGTAGCTAGCTAGCATCGATCGTAGCTAGCATCGTAG"
    "CTAGCATCGAGGCTTAGCGATCGTAGCTAGCATGGCTAGCTAGCTAGCATCGATCGATCG"
)


def load_example():
    st.session_state["primer_input"] = example_seq


col1, col2 = st.columns(2)

with col1:
    st.button("Use example sequence", on_click=load_example)

with col2:
    design = st.button("Design primers", type="primary")


raw_input = st.text_area(
    "Paste your DNA sequence (the target to amplify):",
    height=150,
    placeholder=f"Example: {example_seq[:60]}...",
    key="primer_input"
)


st.subheader("Design Parameters")

col1, col2 = st.columns(2)

with col1:
    primer_len_min, primer_len_max = st.slider(
        "Primer length (bp)",
        min_value=15,
        max_value=30,
        value=(18, 25)
    )

    gc_min, gc_max = st.slider(
        "Target GC content (%)",
        min_value=20,
        max_value=80,
        value=(40, 60)
    )

with col2:
    product_min, product_max = st.slider(
        "Target PCR product size (bp)",
        min_value=50,
        max_value=2000,
        value=(150, 500)
    )

    tm_min, tm_max = st.slider(
        "Target melting temperature (°C)",
        min_value=40,
        max_value=80,
        value=(55, 65)
    )


if design:

    seq = clean_sequence(raw_input)

    if not seq:
        st.error("Please enter a DNA sequence.")
        st.stop()

    invalid = validate_sequence(seq)

    if invalid:
        st.error(
            f"Invalid character(s): {', '.join(invalid)}. "
            "Only A, T, G and C are allowed."
        )
        st.stop()

    if len(seq) < primer_len_min * 2 + 20:
        st.error(
            "Sequence is too short to design a primer pair with these "
            "settings. Enter a longer sequence or reduce the minimum "
            "product size."
        )
        st.stop()

    pairs = design_primer_pairs(
        seq,
        length_range=(primer_len_min, primer_len_max),
        product_range=(product_min, product_max),
        target_gc=(gc_min, gc_max),
        target_tm=(tm_min, tm_max)
    )

    st.session_state["primer_seq"] = seq
    st.session_state["primer_pairs"] = pairs

    if pairs:
        st.success(f"Found {len(pairs)} candidate primer pair(s)!")
    else:
        st.warning(
            "No primer pairs matched these settings. Try widening the "
            "GC%, Tm, or product size ranges."
        )


if "primer_pairs" in st.session_state and st.session_state["primer_pairs"]:

    seq = st.session_state["primer_seq"]
    pairs = st.session_state["primer_pairs"]

    st.divider()

    st.subheader("Candidate Primer Pairs")

    for i, pair in enumerate(pairs, start=1):

        st.markdown(f"### Pair {i} — product size {pair['product_size']} bp")

        col1, col2 = st.columns(2)

        with col1:
            st.markdown("**Forward primer (5' → 3')**")
            st.code(pair["forward_seq"], language="text")
            st.caption(
                f"Position {pair['forward_start']}–{pair['forward_end']} • "
                f"{len(pair['forward_seq'])} bp • "
                f"GC {pair['forward_gc']}% • "
                f"Tm {pair['forward_tm']}°C"
            )

        with col2:
            st.markdown("**Reverse primer (5' → 3')**")
            st.code(pair["reverse_seq"], language="text")
            st.caption(
                f"Position {pair['reverse_start']}–{pair['reverse_end']} • "
                f"{len(pair['reverse_seq'])} bp • "
                f"GC {pair['reverse_gc']}% • "
                f"Tm {pair['reverse_tm']}°C"
            )

        st.caption(f"Tm difference between primers: {pair['tm_diff']}°C")

        st.divider()

    st.subheader("Target Sequence")

    st.code(seq, language="text")

    st.caption(f"Length: {len(seq)} bp")
