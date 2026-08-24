import streamlit as st

from msa_utils import (
    parse_multi_fasta,
    validate_sequences,
    star_alignment,
    consensus_sequence,
    conservation_line,
    pairwise_identity_matrix
)


st.title("Clustal Alignment")

st.write(
    "Paste two or more DNA sequences to produce a multiple sequence "
    "alignment, consensus sequence, and conservation summary."
)

st.caption(
    "This uses a simplified progressive (\"center-star\") alignment "
    "algorithm inspired by ClustalW. It's a lightweight, dependency-free "
    "aligner intended for teaching and quick comparisons, not a drop-in "
    "replacement for the full ClustalW/Clustal Omega tools."
)


example_fasta = """>Seq1
ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCGTGGGGA
>Seq2
ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCATGGGGA
>Seq3
ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCGTGGGA
>Seq4
ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCGTGGGGA
"""


def load_example():
    st.session_state["msa_input"] = example_fasta


col1, col2 = st.columns(2)

with col1:
    st.button("Use example sequences", on_click=load_example)

with col2:
    align = st.button("Align sequences", type="primary")


raw_input = st.text_area(
    "Paste sequences (multi-FASTA, or one sequence per line):",
    height=220,
    placeholder=example_fasta,
    key="msa_input"
)


if align:

    entries = parse_multi_fasta(raw_input)

    if len(entries) < 2:
        st.error("Please enter at least two DNA sequences.")
        st.stop()

    invalid = validate_sequences(entries)

    if invalid:
        details = "; ".join(
            f"{name}: {', '.join(chars)}" for name, chars in invalid.items()
        )
        st.error(
            f"Invalid character(s) found — only A, T, G and C are allowed. "
            f"({details})"
        )
        st.stop()

    names, aligned_seqs = star_alignment(entries)

    st.session_state["msa_names"] = names
    st.session_state["msa_aligned"] = aligned_seqs
    st.session_state["msa_raw_entries"] = entries

    st.success(f"Aligned {len(entries)} sequences successfully!")


if "msa_aligned" in st.session_state:

    names = st.session_state["msa_names"]
    aligned_seqs = st.session_state["msa_aligned"]
    raw_entries = st.session_state["msa_raw_entries"]

    consensus = consensus_sequence(aligned_seqs)
    conservation = conservation_line(aligned_seqs)

    identical_cols = conservation.count("*")
    alignment_length = len(consensus)
    overall_conservation = (
        round(identical_cols / alignment_length * 100, 2)
        if alignment_length else 0.0
    )

    st.divider()

    st.subheader("Alignment Summary")

    col1, col2, col3 = st.columns(3)

    with col1:
        st.metric("Sequences", len(names))

    with col2:
        st.metric("Alignment length", f"{alignment_length} bp")

    with col3:
        st.metric("Fully conserved columns", f"{overall_conservation}%")

    st.divider()

    st.subheader("Multiple Sequence Alignment")

    st.caption("'*' marks columns where every sequence agrees, with no gaps.")

    name_width = max(len(name) for name in names + ["Consensus"]) + 2

    chunk_size = 60

    for start in range(0, alignment_length, chunk_size):
        end = start + chunk_size

        lines = []

        for name, seq in zip(names, aligned_seqs):
            lines.append(f"{name.ljust(name_width)}{seq[start:end]}")

        lines.append(f"{''.ljust(name_width)}{conservation[start:end]}")
        lines.append(f"{'Consensus'.ljust(name_width)}{consensus[start:end]}")

        st.code(
            f"Position {start + 1} - {min(end, alignment_length)}\n\n"
            + "\n".join(lines),
            language="text"
        )

    st.divider()

    st.subheader("Consensus Sequence")

    st.code(consensus, language="text")

    st.divider()

    st.subheader("Pairwise Identity Matrix")

    matrix = pairwise_identity_matrix(names, aligned_seqs)

    header_row = "| " + " | ".join([" "] + names) + " |"
    separator_row = "|" + "---|" * (len(names) + 1)

    table_lines = [header_row, separator_row]

    for name, row in zip(names, matrix):
        cells = " | ".join(f"{value:.1f}%" for value in row)
        table_lines.append(f"| **{name}** | {cells} |")

    st.markdown("\n".join(table_lines))

    st.divider()

    st.subheader("Original Sequences")

    for name, seq in raw_entries:
        st.markdown(f"**{name}** ({len(seq)} bp)")
        st.code(seq, language="text")
