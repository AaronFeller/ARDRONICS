# JOSS submission notes

This is an early manuscript draft, not yet a submission-ready claim of biological validation.

## Files

- `paper.md`: JOSS-format source with YAML metadata. Current JOSS documentation requires Markdown.
- `paper.tex`: LaTeX companion requested for sharing and editing.
- `paper.bib`: shared bibliography.
- `figures/workflow.png` and `workflow.pdf`: original workflow figure.
- `make_figure.py`: reproducibly regenerates the figure.

## Before submission

1. Add ORCID identifiers as appropriate. Affiliations are currently UT Austin Interdisciplinary Life Sciences Graduate Programs for Aaron L. Feller and ASU School of Life Sciences for Hinsby Cadillo-Quiroz.
2. Add the public source-repository URL and a versioned Zenodo (or equivalent) archive DOI.
3. Review authorship, conflicts of interest, and the AI-usage disclosure for accuracy.
4. Conduct the planned blinded biological validation and add actual results, or expect editors to request stronger evidence of research impact.
5. Confirm that any figures or example datasets made public have redistribution permission. The workflow figure is original; no third-party gel image is used in the manuscript.
6. Review every reference and compile the paper through the official Open Journals GitHub Action or `openjournals/inara` Docker image.
7. Tag the software version reviewed by JOSS and ensure the archived repository includes tests, documentation, license, and manuscript.

## Current verification snapshot

- 29 tests passing
- 89% statement coverage overall; 92% digestion and 97% matching
- Ruff static checks passing
- Wheel and source distribution builds passing
- Installed dependency compatibility check passing

These are software-verification results, not estimates of organism-identification accuracy.
