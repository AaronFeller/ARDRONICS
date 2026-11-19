.PHONY: all latex markdown figure clean

all: latex markdown

latex: figure
	tectonic paper.tex

markdown: figure
	pandoc paper.md --citeproc --pdf-engine=tectonic --variable=author:"Aaron L. Feller and Hinsby Cadillo-Quiroz" -o paper_markdown.pdf

figure:
	../.venv/bin/python make_figure.py

clean:
	rm -f paper.pdf paper_markdown.pdf
