# 🧬 DNA Sequence Analysis Toolkit

A single-file Streamlit app. Paste a DNA sequence and get:

1. **Length** of the sequence
2. **Complementary strand** (and reverse complement)
3. **A/T/G/C composition** (percentages + GC content)
4. **RNA sequence** (T → U transcription)
5. **Protein sequence** (translation, reading frame 1)
6. **Predicted 3D protein structure** (folded via the free ESMFold API, viewed in-browser)
7. **Protein BLAST** (searched live against NCBI's `nr` database via `blastp`)

## Setup

```bash
python -m venv venv
source venv/bin/activate      # on Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## Run

```bash
streamlit run app.py
```

This opens the app at `http://localhost:8501`.

## Notes & limitations

- **Internet required** for steps 6 and 7 — they call external, free public
  APIs (ESMFold / ESM Metagenomic Atlas, and NCBI BLAST). Steps 1-5 work
  fully offline.
- **ESMFold structure prediction** works best on sequences under ~400 amino
  acids (a limit of the free public API). It folds only up to the first
  stop codon (`*`).
- **NCBI BLAST** can take 1-5 minutes per search because it runs on NCBI's
  servers and is queued. Avoid firing off many requests back-to-back — NCBI
  rate-limits automated queries. For heavy/production use, consider running
  a local BLAST installation against a local database instead.
- Input is cleaned automatically: FASTA header lines (starting with `>`),
  whitespace, and line breaks are stripped. Any character outside `A T G C`
  will raise a validation error.
- Translation reads only frame 1, starting at the first base — it does not
  search for the first `ATG` start codon. If your sequence includes 5'
  untranslated region before the start codon, trim it first for a cleaner
  protein output.

## Deploying it as a real website

To make this reachable at a URL rather than just `localhost`, the easiest
options are:

- **Streamlit Community Cloud** (free) — push this folder to a GitHub repo
  and connect it at share.streamlit.io.
- Any VM/server: run `streamlit run app.py --server.port 80` behind a
  reverse proxy (e.g. nginx), or containerize it with Docker.

## Extending it

- Swap ESMFold for a local ColabFold/AlphaFold install if you need longer
  sequences or offline structure prediction.
- Swap NCBI's remote BLAST for a local BLAST+ install (`blastp` binary +
  a downloaded database) if you need speed or to avoid rate limits.
- Add multi-frame translation (frames 1-3, plus reverse strand) if you want
  full 6-frame translation like ExPASy Translate.
