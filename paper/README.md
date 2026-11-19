# Manuscript builds

`paper.md` is the current JOSS-format source. `paper.tex` is the requested LaTeX companion for sharing and editing; keep both synchronized when revising content.

Build locally with Pandoc and Tectonic installed:

```bash
cd paper
make latex      # paper.pdf with BibTeX references
make markdown   # paper_markdown.pdf using the JOSS Markdown source
make figure     # regenerate the original workflow figure
```

For actual JOSS submission, also compile `paper.md` using the official Open Journals GitHub Action or `openjournals/inara` Docker image. See `SUBMISSION_NOTES.md` for remaining author and validation tasks.
