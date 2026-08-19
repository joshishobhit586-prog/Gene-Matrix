import streamlit as st


st.set_page_config(
    page_title="DNA Analysis Toolkit",
    page_icon="🧬",
    layout="wide"
)


st.title("🧬 DNA Sequence Analysis Toolkit")

st.markdown(
    """
    Welcome to the **DNA Sequence Analysis Toolkit**.

    This application allows you to analyze a DNA sequence and explore:

    - 🧬 Basic sequence information
    - 🔗 Complementary and reverse-complementary strands
    - 📊 A/T/G/C composition
    - 🧪 RNA transcription
    - 🧬 Protein translation
    - 🧫 Predicted 3D protein structure

    ### How to use

    1. Go to **Sequence Analysis** in the sidebar.
    2. Enter your DNA sequence.
    3. Click **Analyze Sequence**.
    4. Use the sidebar to explore the different analyses.
    """
)


st.divider()

st.subheader("Example DNA sequence")

example_seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

st.code(example_seq, language="text")

st.info(
    "Start by opening **🧬 Sequence Analysis** from the sidebar."
)