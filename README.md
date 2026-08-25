# Gene Matrix

Gene Matrix is an interactive bioinformatics toolkit built with Python and Streamlit. It provides tools for analyzing DNA sequences, exploring nucleotide composition, performing transcription and translation, comparing sequences, creating multiple sequence alignments, designing PCR primers, and predicting protein structures.


## Live Application
Try Gene Matrix here:
[Open Gene Matrix](https://dna-sequenza.streamlit.app?utm_source=chatgpt.com)


## Features

### Sequence Analysis

Analyze DNA sequences to obtain:

* Sequence length
* Complementary DNA strand
* Reverse complement
* Sequence validation and cleaning

### Base Composition

Calculate:

* Adenine (A)
* Thymine (T)
* Guanine (G)
* Cytosine (C)
* GC content

### Transcription and Translation

Convert DNA into RNA and translate the resulting sequence into its corresponding protein sequence.

### Protein Structure Prediction

Predict protein structures using the ESMFold API.

Features include:

* Protein sequence input
* Structure prediction
* Predicted structure visualization
* API error handling

### DNA Sequence Comparison

Compare two DNA sequences and analyze:

* Sequence lengths
* Matches
* Mismatches
* Gaps
* Sequence similarity
* Alignment results

### Multiple Sequence Alignment

Perform multiple sequence alignment and generate:

* Aligned sequences
* Consensus sequence
* Conserved positions
* Pairwise identity matrix
* Alignment statistics

### PCR Primer Designer

Generate and evaluate candidate forward and reverse primers based on:

* Primer length
* GC content
* Melting temperature
* Primer compatibility
* Tm difference
* Primer ranking

## Technologies Used

* Python
* Streamlit
* Biopython
* Requests
* ESMFold API

## Project Structure

```text
gene-matrix/
│
├── Home.py
├── dna_utils.py
├── msa_utils.py
├── primer_utils.py
├── requirements.txt
│
├── pages/
│   ├── 1_Sequence_Analysis.py
│   ├── 2_Composition.py
│   ├── 3_Translation.py
│   ├── 4_Protein_Structure.py
│   ├── 5_Sequence_Comparison.py
│   ├── 6_Clustal_Alignment.py
│   └── 7_PCR_Primer_Designer.py
│
└── README.md
```

## Installation

Clone the repository:

```bash
git clone https://github.com/joshishobhit586-prog/Gene-Matrix.git
```

Navigate to the project directory:

```bash
cd Gene-Matrix
```

Create a virtual environment:

### Linux/macOS

```bash
python -m venv venv
source venv/bin/activate
```

### Windows

```bash
python -m venv venv
venv\Scripts\activate
```

Install the dependencies:

```bash
pip install -r requirements.txt
```

## Running the Application

Run the application with:

```bash
streamlit run Home.py
```

## Project Architecture

### dna_utils.py

Handles core DNA operations including:

* Sequence cleaning
* Sequence validation
* Complement generation
* Reverse complement generation
* Sequence comparison

### msa_utils.py

Handles:

* Multiple sequence alignment
* Pairwise alignment
* Consensus sequence generation
* Identity calculations

### primer_utils.py

Handles:

* Primer generation
* GC content calculation
* Melting temperature estimation
* Primer pair ranking

### pages/

Contains the individual Streamlit pages for each bioinformatics tool.

## Disclaimer

Gene Matrix is designed for educational and exploratory purposes. The analyses, alignments, primer suggestions, and structure predictions should not be used as a substitute for professional laboratory validation or clinical analysis.
