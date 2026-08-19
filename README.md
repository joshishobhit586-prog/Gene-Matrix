# 🧬 DNA Sequence Analysis Toolkit

A single-file Streamlit app. Paste a DNA sequence and get:

1. **Length** of the sequence
2. **Complementary strand** (and reverse complement)
3. **A/T/G/C composition** (percentages + GC content)
4. **RNA sequence** (T → U transcription)
5. **Protein sequence** (translation, reading frame 1)
6. **Predicted 3D protein structure** (folded via the free ESMFold API, viewed in-browser)

## Setup

```bash
python -m venv venv
source venv/bin/activate      # on Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## Run

```url
https://gene-matrix-uuwjruximt3u6bn67kuqvp.streamlit.app/Sequence_Analysis```

This opens the Website

## Notes & limitations

- **Internet required** for steps 6  — it calls an external, free public
  API (ESMFold / ESM Metagenomic Atlas). Steps 1-5 work
  fully offline.
- **ESMFold structure prediction** works best on sequences under ~400 amino
  acids (a limit of the free public API). It folds only up to the first
  stop codon (`*`).
- Input is cleaned automatically: FASTA header lines (starting with `>`),
  whitespace, and line breaks are stripped. Any character outside `A T G C`
  will raise a validation error.
- Translation reads only frame 1, starting at the first base — it does not
  search for the first `ATG` start codon. If your sequence includes 5'
  untranslated region before the start codon, trim it first for a cleaner
  protein output.

## Deploying it as a real website

To make this reachable at a URL rather than just `localhost`, the easiest
option is:

- **Streamlit Community Cloud** (free) — push this folder to a GitHub repo
  and connect it at share.streamlit.io.

