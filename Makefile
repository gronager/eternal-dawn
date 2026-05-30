MAIN := book
BUILD := build
SIM := sims
VENV := $(SIM)/.venv
PYTHON := $(VENV)/bin/python
FIGSCRIPTS := $(wildcard figures/scripts/ch*.py)

.PHONY: all pdf figures epub sim-install sim-test clean distclean

all: pdf

pdf:
	@mkdir -p $(BUILD)
	latexmk -pdf $(MAIN).tex

# Regenerate every chapter figure (PDFs land in figures/pdf/).
figures: sim-install
	@for f in $(FIGSCRIPTS); do echo ">> $$f"; $(PYTHON) $$f || exit 1; done

# EPUB3 via LaTeXML. latexmlc converts math to MathML and rasterises the
# \includegraphics PDFs (needs Ghostscript/ImageMagick on PATH). Install LaTeXML
# with e.g. 'cpanm LaTeXML' or your distro's 'latexml' package.
epub:
	@command -v latexmlc >/dev/null 2>&1 || { \
		echo "latexmlc not found. Install LaTeXML (e.g. 'cpanm LaTeXML' or apt/brew 'latexml'),"; \
		echo "plus Ghostscript/ImageMagick for figure conversion."; exit 1; }
	@mkdir -p $(BUILD)
	latexmlc $(MAIN).tex --dest=$(BUILD)/$(MAIN).epub \
		--splitat=chapter --timeout=600

sim-install: $(VENV)
$(VENV):
	cd $(SIM) && uv venv && uv pip install -e ".[dev]"

sim-test: sim-install
	cd $(SIM) && .venv/bin/pytest -q

clean:
	latexmk -c

distclean:
	latexmk -C
	rm -rf $(BUILD) $(VENV)
