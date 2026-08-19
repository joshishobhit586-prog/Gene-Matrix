import streamlit as st

from dna_utils import (
    clean_sequence,
    validate_sequence,
    complementary_strand,
    reverse_complement
)


st.title("🧬 Sequence Analysis")

st.write(
    "Enter a DNA sequence to calculate its length, "
    "complementary strand, and reverse complement."
)


example_seq = "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCGTGGGGA"


raw_input = st.text_area(
    "Paste your DNA sequence:",
    height=150,
    placeholder=f"Example: {example_seq}",
    key="dna_input"
)


col1, col2 = st.columns(2)

with col1:
    use_example = st.button("Use example sequence")

with col2:
    analyze = st.button(
        "🔍 Analyze sequence",
        type="primary"
    )


def load_example():
    st.session_state["dna_input"] = example_seq

st.button(
    "Use example sequence",
    on_click=load_example
)

'''
old and not wroking
if use_example:
    st.session_state["dna_input"] = example_seq
    st.rerun()
'''

if analyze:

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

    # Store results in session state
    st.session_state["seq"] = seq

    st.session_state["comp"] = complementary_strand(seq)

    st.session_state["rev_comp"] = reverse_complement(seq)

    st.session_state["composition"] = {
        base: round(seq.count(base) / len(seq) * 100, 2)
        for base in "ATGC"
    }

    st.session_state["rna"] = seq.replace("T", "U")

    from dna_utils import translate_dna

    st.session_state["protein"] = translate_dna(seq)

    st.success("Sequence analyzed successfully!")


# Display results if a sequence has been analyzed

if "seq" in st.session_state:

    seq = st.session_state["seq"]

    st.divider()

    st.subheader("Basic Information")

    st.metric(
        "Sequence length",
        f"{len(seq)} bp"
    )

    st.markdown("**DNA sequence:**")

    st.code(
        seq,
        language="text"
    )

    st.divider()

    st.subheader("Complementary Strand")

    col1, col2 = st.columns(2)

    with col1:

        st.markdown("**5' → 3' Original**")

        st.code(
            seq,
            language="text"
        )

        st.markdown("**3' → 5' Complement**")

        st.code(
            st.session_state["comp"],
            language="text"
        )

    with col2:

        st.markdown("**5' → 3' Reverse Complement**")

        st.code(
            st.session_state["rev_comp"],
            language="text"
        )