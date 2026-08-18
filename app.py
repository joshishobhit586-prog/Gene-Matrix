"""
DNA Sequence Analysis Toolkit

A Streamlit web app that takes a DNA sequence and returns:
  1. Sequence length
  2. Complementary strand (and reverse complement)
  3. A/T/G/C composition (percentages)
  4. RNA sequence (transcription)
  5. Protein sequence (translation)
  6. Predicted 3D protein structure (via the ESMFold API)
"""

import re
import requests
import streamlit as st
from Bio.Seq import Seq


# Page setup

st.set_page_config(page_title="DNA Analysis Toolkit", page_icon="🧬", layout="wide")
st.title("🧬 DNA Sequence Analysis Toolkit")
st.caption(
    "Enter a DNA sequence to get its length, complementary strand, base "
    "composition, RNA transcript, protein translation, a predicted 3D "
    "protein structure"
)


# Helpers


VALID_BASES = set("ATGC")
COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


def clean_sequence(raw: str) -> str:
    """Uppercase and strip whitespace/newlines/numbers (e.g. pasted FASTA)."""
    seq = raw.strip().upper()
    lines = seq.splitlines()
    # drop a FASTA header line if present
    lines = [ln for ln in lines if not ln.startswith(">")]
    seq = "".join(lines)
    seq = re.sub(r"[^A-Z]", "", seq)  # remove spaces/digits/punctuation
    return seq


def validate_sequence(seq: str):
    invalid = sorted(set(seq) - VALID_BASES)
    return invalid


def complementary_strand(seq: str) -> str:
    return "".join(COMPLEMENT[b] for b in seq)


def reverse_complement(seq: str) -> str:
    return complementary_strand(seq)[::-1]


def base_composition(seq: str) -> dict:
    n = len(seq)
    return {b: round(seq.count(b) / n * 100, 2) for b in "ATGC"} if n else {}


def to_rna(seq: str) -> str:
    return seq.replace("T", "U")


def translate_dna(seq: str) -> str:
    """Translate DNA -> protein using Biopython (frame 1)."""
    coding_len = len(seq) - (len(seq) % 3)
    protein = str(Seq(seq[:coding_len]).translate(to_stop=False))
    return protein

#Input

example_seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

raw_input = st.text_area(
    "Paste a DNA sequence (FASTA header lines starting with '>' are ignored):",
    height=150,
    placeholder=f"e.g. {example_seq}",
)

col_a, col_b = st.columns([1, 5])
with col_a:
    use_example = st.button("Use example sequence")
if use_example:
    raw_input = example_seq
    st.session_state["raw_input"] = example_seq

analyze = st.button("🔍 Analyze sequence", type="primary")

if analyze:
    seq = clean_sequence(raw_input)

    if not seq:
        st.error("Please enter a DNA sequence.")
        st.stop()

    invalid = validate_sequence(seq)
    if invalid:
        st.error(
            f"Sequence contains invalid character(s): {', '.join(invalid)}. "
            "Only A, T, G, C are allowed."
        )
        st.stop()

    st.session_state["seq"] = seq
    st.session_state["comp"] = complementary_strand(seq)
    st.session_state["rev_comp"] = reverse_complement(seq)
    st.session_state["composition"] = base_composition(seq)
    st.session_state["rna"] = to_rna(seq)
    st.session_state["protein"] = translate_dna(seq)


#Output

if "seq" in st.session_state:
    seq = st.session_state["seq"]

    st.divider()
    st.subheader("1. Basic Info")
    st.metric("Sequence length (bp)", len(seq))
    st.code(seq, language="text")

    st.subheader("2. Complementary Strand")
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("**5'→3' (original)**")
        st.code(seq, language="text")
        st.markdown("**3'→5' (complement)**")
        st.code(st.session_state["comp"], language="text")
    with c2:
        st.markdown("**5'→3' (reverse complement)**")
        st.code(st.session_state["rev_comp"], language="text")

    st.subheader("3. Base Composition (A/T/G/C %)")
    comp_pct = st.session_state["composition"]
    b1, b2, b3, b4 = st.columns(4)
    for col, base in zip((b1, b2, b3, b4), "ATGC"):
        col.metric(base, f"{comp_pct[base]}%")
    st.bar_chart(comp_pct)
    gc_content = comp_pct["G"] + comp_pct["C"]
    st.caption(f"GC content: {gc_content:.2f}%  |  AT content: {100 - gc_content:.2f}%")

    st.subheader("4. RNA Sequence (Transcription)")
    st.code(st.session_state["rna"], language="text")

    st.subheader("5. Protein Sequence (Translation, frame 1)")
    protein = st.session_state["protein"]
    st.code(protein, language="text")
    if len(seq) % 3 != 0:
        st.caption(
            f"Note: sequence length isn't a multiple of 3 — the last "
            f"{len(seq) % 3} base(s) were dropped for translation."
        )
    st.caption("'*' marks a stop codon.")

    clean_protein = protein.replace("*", "")
#-------------------------------------------------------
    st.divider()
    st.subheader("6. Predicted 3D Protein Structure")
    st.caption(
        "Uses the free ESMFold API (ESM Metagenomic Atlas) to predict a 3D "
        "structure from the protein sequence. Best for sequences under ~400 "
        "residues; requires no stop codons ('*') in the sequence sent."
    )

    if not clean_protein:
        st.warning("No protein sequence to fold (empty after removing stop codons).")
    elif "*" in protein[:-1]:
        st.warning(
            "This protein contains an internal stop codon, so only the "
            "sequence up to the first stop will be folded."
        )

    fold_seq = protein.split("*")[0]  # sequence up to first stop codon

    if st.button("🧫 Predict 3D structure with ESMFold"):
        if len(fold_seq) < 2:
            st.error("Protein sequence too short to fold.")
        elif len(fold_seq) > 400:
            st.error(
                f"Sequence is {len(fold_seq)} residues — the free ESMFold API "
                "caps at 400. Try a shorter sequence."
            )
        else:
            with st.spinner("Contacting ESMFold API and predicting structure..."):
                try:
                    resp = requests.post(
                        "https://api.esmatlas.com/foldSequence/v1/pdb/",
                        data=fold_seq,
                        timeout=180,
                    )
                    if resp.status_code == 200 and resp.text.strip().startswith(
                        ("HEADER", "ATOM", "REMARK")
                    ):
                        pdb_text = resp.text
                        st.session_state["pdb_text"] = pdb_text
                        st.success("Structure predicted.")
                    else:
                        st.error(
                            f"ESMFold API returned an error (status "
                            f"{resp.status_code}). Try again later."
                        )
                except requests.RequestException as e:
                    st.error(f"Could not reach ESMFold API: {e}")

    if "pdb_text" in st.session_state:
        try:
            import py3Dmol
            from stmol import showmol

            view = py3Dmol.view(width=800, height=500)
            view.addModel(st.session_state["pdb_text"], "pdb")
            view.setStyle({"cartoon": {"color": "spectrum"}})
            view.zoomTo()
            showmol(view, height=500, width=800)
        except ImportError:
            st.info(
                "Install `py3Dmol` and `stmol` to view the structure in-app "
                "(`pip install py3Dmol stmol`). You can still download the "
                "raw PDB file below."
            )
        st.download_button(
            "⬇️ Download predicted structure (.pdb)",
            data=st.session_state["pdb_text"],
            file_name="predicted_structure.pdb",
            mime="chemical/x-pdb",
        )
else:
    st.info("Enter a DNA sequence above and click **Analyze sequence** to begin.")
