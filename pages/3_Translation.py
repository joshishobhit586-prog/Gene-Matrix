import streamlit as st


st.title("🧪 Transcription & Translation")


if "seq" not in st.session_state:

    st.warning(
        "No DNA sequence has been analyzed yet."
    )

    st.info(
        "Go to **🧬 Sequence Analysis** and analyze a sequence first."
    )

    st.stop()


seq = st.session_state["seq"]

rna = st.session_state["rna"]

protein = st.session_state["protein"]


# RNA

st.subheader("1. RNA Sequence")

st.markdown(
    "**DNA → RNA transcription**"
)

st.code(
    rna,
    language="text"
)


st.divider()


# Protein

st.subheader("2. Protein Sequence")

st.markdown(
    "**DNA → Protein translation (Frame 1)**"
)

st.code(
    protein,
    language="text"
)


if len(seq) % 3 != 0:

    remainder = len(seq) % 3

    st.info(
        f"The DNA sequence is not a multiple of 3. "
        f"The final {remainder} base(s) were not translated."
    )


st.caption(
    "'*' represents a stop codon."
)


st.divider()


st.subheader("Translation Information")


col1, col2, col3 = st.columns(3)


with col1:

    st.metric(
        "DNA length",
        f"{len(seq)} bp"
    )


with col2:

    st.metric(
        "RNA length",
        f"{len(rna)} bases"
    )


with col3:

    st.metric(
        "Protein length",
        f"{len(protein)} residues"
    )