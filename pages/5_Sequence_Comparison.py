import streamlit as st

from dna_utils import (
    clean_sequence,
    validate_sequence,
    compare_sequences
)


st.title("🔬 DNA Sequence Comparison")

st.write(
    "Compare two DNA sequences using global sequence alignment."
)


# --------------------------------------------------
# Example sequences
# --------------------------------------------------

example_seq_1 = (
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCT"
    "GCCGTTACTGCCCGTGGGGA"
)

example_seq_2 = (
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCT"
    "GCCGTTACTGCCCATGGGGA"
)


# --------------------------------------------------
# Load example sequences
# --------------------------------------------------

def load_examples():

    st.session_state["comparison_seq1"] = (
        example_seq_1
    )

    st.session_state["comparison_seq2"] = (
        example_seq_2
    )


# --------------------------------------------------
# Input section
# --------------------------------------------------

st.subheader("Enter Two DNA Sequences")

col1, col2 = st.columns(2)

with col1:

    st.text_area(
        "🧬 DNA Sequence 1",
        height=200,
        placeholder="Enter first DNA sequence...",
        key="comparison_seq1"
    )


with col2:

    st.text_area(
        "🧬 DNA Sequence 2",
        height=200,
        placeholder="Enter second DNA sequence...",
        key="comparison_seq2"
    )


# --------------------------------------------------
# Buttons
# --------------------------------------------------

button_col1, button_col2 = st.columns(2)

with button_col1:

    st.button(
        "Use example sequences",
        on_click=load_examples
    )


with button_col2:

    compare = st.button(
        "🔬 Compare sequences",
        type="primary",
        use_container_width=True
    )


# --------------------------------------------------
# Compare sequences
# --------------------------------------------------

if compare:

    raw_seq1 = (
        st.session_state["comparison_seq1"]
    )

    raw_seq2 = (
        st.session_state["comparison_seq2"]
    )

    seq1 = clean_sequence(raw_seq1)

    seq2 = clean_sequence(raw_seq2)


    # Check for empty sequences

    if not seq1:

        st.error(
            "Please enter DNA Sequence 1."
        )

        st.stop()


    if not seq2:

        st.error(
            "Please enter DNA Sequence 2."
        )

        st.stop()


    # Validate sequence 1

    invalid1 = validate_sequence(seq1)

    if invalid1:

        st.error(
            "Sequence 1 contains invalid character(s): "
            f"{', '.join(invalid1)}. "
            "Only A, T, G and C are allowed."
        )

        st.stop()


    # Validate sequence 2

    invalid2 = validate_sequence(seq2)

    if invalid2:

        st.error(
            "Sequence 2 contains invalid character(s): "
            f"{', '.join(invalid2)}. "
            "Only A, T, G and C are allowed."
        )

        st.stop()


    # Perform comparison

    results = compare_sequences(
        seq1,
        seq2
    )


    # Store results

    st.session_state["comparison_results"] = (
        results
    )


    st.success(
        "Sequences compared successfully!"
    )


# --------------------------------------------------
# Display results
# --------------------------------------------------

if "comparison_results" in st.session_state:

    results = (
        st.session_state["comparison_results"]
    )


    st.divider()

    st.header("📊 Comparison Results")


    # --------------------------------------------------
    # Main metrics
    # --------------------------------------------------

    col1, col2, col3, col4 = st.columns(4)


    with col1:

        st.metric(
            "Sequence Identity",
            f"{results['identity']}%"
        )


    with col2:

        st.metric(
            "Matches",
            results["matches"]
        )


    with col3:

        st.metric(
            "Mismatches",
            results["mismatches"]
        )


    with col4:

        st.metric(
            "Gaps",
            results["gaps"]
        )


    # --------------------------------------------------
    # Sequence lengths
    # --------------------------------------------------

    st.divider()

    st.subheader("📏 Sequence Lengths")


    col1, col2, col3 = st.columns(3)


    with col1:

        st.metric(
            "Sequence 1",
            f"{results['length_1']} bp"
        )


    with col2:

        st.metric(
            "Sequence 2",
            f"{results['length_2']} bp"
        )


    with col3:

        st.metric(
            "Aligned Length",
            f"{results['aligned_length']} bp"
        )


    # --------------------------------------------------
    # Visual alignment
    # --------------------------------------------------

    st.divider()

    st.subheader("🧬 Sequence Alignment")

    st.caption(
        "| = Match   •   * = Mismatch   •   "
        "space = Gap"
    )


    aligned_seq1 = (
        results["aligned_seq1"]
    )

    aligned_seq2 = (
        results["aligned_seq2"]
    )

    comparison_line = (
        results["comparison_line"]
    )


    chunk_size = 60


    for start in range(
        0,
        len(aligned_seq1),
        chunk_size
    ):

        end = start + chunk_size


        seq1_chunk = (
            aligned_seq1[start:end]
        )

        match_chunk = (
            comparison_line[start:end]
        )

        seq2_chunk = (
            aligned_seq2[start:end]
        )


        st.code(
            f"""Position {start + 1} - {min(end, len(aligned_seq1))}

Sequence 1: {seq1_chunk}
            {match_chunk}
Sequence 2: {seq2_chunk}
""",
            language="text"
        )


    # --------------------------------------------------
    # Original sequences
    # --------------------------------------------------

    st.divider()

    st.subheader("Original Sequences")


    col1, col2 = st.columns(2)


    with col1:

        st.markdown(
            "**DNA Sequence 1**"
        )

        st.code(
            results["seq1"],
            language="text"
        )


    with col2:

        st.markdown(
            "**DNA Sequence 2**"
        )

        st.code(
            results["seq2"],
            language="text"
        )


    # --------------------------------------------------
    # Interpretation
    # --------------------------------------------------

    st.divider()

    st.subheader("🧪 Interpretation")

    identity = results["identity"]


    if identity == 100:

        st.success(
            "The two DNA sequences are identical."
        )


    elif identity >= 95:

        st.success(
            "The sequences are highly similar."
        )


    elif identity >= 70:

        st.info(
            "The sequences have moderate similarity."
        )


    else:

        st.warning(
            "The sequences have low similarity."
        )