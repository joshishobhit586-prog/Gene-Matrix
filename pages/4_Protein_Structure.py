import time
import requests
import streamlit as st


st.title("🧫 3D Protein Structure")


if "protein" not in st.session_state:

    st.warning(
        "No protein sequence has been generated yet."
    )

    st.info(
        "Go to **🧬 Sequence Analysis** and analyze a DNA sequence first."
    )

    st.stop()


protein = st.session_state["protein"]

clean_protein = protein.replace("*", "")


st.subheader("Protein Sequence")

st.code(
    protein,
    language="text"
)


st.caption(
    "ESMFold works with the protein sequence generated from your DNA."
)


if not clean_protein:

    st.warning(
        "No protein sequence is available for structure prediction."
    )

    st.stop()


if "*" in protein[:-1]:

    st.warning(
        "The protein contains an internal stop codon. "
        "Only the sequence before the first stop codon will be folded."
    )


fold_seq = protein.split("*")[0]


st.write(
    f"Protein length used for folding: **{len(fold_seq)} residues**"
)


if len(fold_seq) < 2:

    st.error(
        "Protein sequence is too short to fold."
    )

    st.stop()


if len(fold_seq) > 400:

    st.error(
        f"Protein sequence is {len(fold_seq)} residues long. "
        "The free ESMFold API is limited to approximately 400 residues."
    )

    st.stop()


st.divider()


predict = st.button(
    "🧫 Predict 3D Structure with ESMFold",
    type="primary"
)


if predict:

    RETRYABLE_STATUSES = {
        502,
        503,
        504
    }

    MAX_ATTEMPTS = 5

    pdb_text = None

    last_status = None

    last_error = None


    status_box = st.empty()


    for attempt in range(1, MAX_ATTEMPTS + 1):

        status_box.info(
            f"Contacting ESMFold API "
            f"(attempt {attempt}/{MAX_ATTEMPTS})..."
        )

        try:

            response = requests.post(
                "https://api.esmatlas.com/foldSequence/v1/pdb/",
                data=fold_seq,
                timeout=180
            )


            if (
                response.status_code == 200
                and response.text.strip().startswith(
                    ("HEADER", "ATOM", "REMARK")
                )
            ):

                pdb_text = response.text

                break


            last_status = response.status_code


            if response.status_code not in RETRYABLE_STATUSES:

                break


        except requests.RequestException as error:

            last_error = error


        if attempt < MAX_ATTEMPTS:

            time.sleep(5 * attempt)


    status_box.empty()


    if pdb_text:

        st.session_state["pdb_text"] = pdb_text

        st.success(
            "3D structure predicted successfully!"
        )

    else:

        st.error(
            "The ESMFold service did not return a structure."
        )

        if last_status:

            st.caption(
                f"API status code: {last_status}"
            )

        if last_error:

            st.caption(
                f"Connection error: {last_error}"
            )


# Display structure if it exists

if "pdb_text" in st.session_state:

    st.divider()

    st.subheader("Predicted 3D Structure")


    try:

        import py3Dmol
        from stmol import showmol


        view = py3Dmol.view(
            width=800,
            height=500
        )


        view.addModel(
            st.session_state["pdb_text"],
            "pdb"
        )


        view.setStyle(
            {
                "cartoon": {
                    "color": "spectrum"
                }
            }
        )


        view.zoomTo()


        showmol(
            view,
            height=500,
            width=800
        )


    except ImportError:

        st.info(
            "Install py3Dmol and stmol to display "
            "the structure inside the application."
        )


    st.download_button(
        "⬇️ Download predicted structure (.pdb)",
        data=st.session_state["pdb_text"],
        file_name="predicted_structure.pdb",
        mime="chemical/x-pdb"
    )