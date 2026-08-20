import streamlit as st


st.title("Base Composition")


if "seq" not in st.session_state:

    st.warning(
        "No DNA sequence has been analyzed yet."
    )

    st.info(
        "Go to **Sequence Analysis** and analyze a sequence first."
    )

    st.stop()


seq = st.session_state["seq"]

composition = st.session_state["composition"]


st.subheader("A / T / G / C Composition")


col1, col2, col3, col4 = st.columns(4)


with col1:
    st.metric(
        "Adenine (A)",
        f"{composition['A']}%"
    )


with col2:
    st.metric(
        "Thymine (T)",
        f"{composition['T']}%"
    )


with col3:
    st.metric(
        "Guanine (G)",
        f"{composition['G']}%"
    )


with col4:
    st.metric(
        "Cytosine (C)",
        f"{composition['C']}%"
    )


st.divider()


st.subheader("Base Composition Chart")

st.bar_chart(composition)


gc_content = composition["G"] + composition["C"]

at_content = composition["A"] + composition["T"]


st.divider()


col1, col2 = st.columns(2)


with col1:

    st.metric(
        "GC Content",
        f"{gc_content:.2f}%"
    )


with col2:

    st.metric(
        "AT Content",
        f"{at_content:.2f}%"
    )


st.divider()


st.subheader("Analyzed Sequence")

st.code(
    seq,
    language="text"
)
