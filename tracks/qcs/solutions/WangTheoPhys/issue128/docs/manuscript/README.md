# Long-form Issue 128 manuscript

This directory contains the full English research-paper treatment of the
machine-certified Issue 128 result. The authoritative current numbers remain
in `../../certificates/issue128-d5-integrated-certificate.json` and its hashed
D4/D5 sidecars; the PDF is a human-readable derived artifact. The original
97-step certificate remains frozen as a separately auditable predecessor.

## Build

From this directory:

```bash
make pdf
```

The final file is:

```text
output/pdf/issue128-proof-carrying-trotter-paper.pdf
```

The build first validates headline claims against the frozen JSON, regenerates
the four vector PDF figures, and then runs `latexmk` with BibTeX until all
cross-references settle.

Equivalent commands are:

```bash
python3 scripts/validate_claims.py
python3 scripts/generate_figures.py
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  -jobname=issue128-proof-carrying-trotter-paper \
  -outdir=output/pdf main.tex
```

## Verification levels

From the Issue 128 root:

```bash
pytest -q
PYTHONPATH=src python3 scripts/verify.py \
  certificates/issue128-d5-integrated-certificate.json
python3 scripts/reference_verify.py \
  certificates/issue128-d5-integrated-certificate.json
```

The expensive independent replay is optional:

```bash
PYTHONPATH=src python3 scripts/verify.py \
  certificates/issue128-d5-integrated-certificate.json --deep
```

See `CITATION_AUDIT.md` for the reference-validation result. Generated LaTeX
auxiliary files are ignored; the named final PDF is retained.
