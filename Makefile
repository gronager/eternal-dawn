MAIN := book
BUILD := build
SIM := sims
VENV := $(SIM)/.venv
PYTHON := $(VENV)/bin/python
FIGSCRIPTS := $(wildcard figures/scripts/ch*.py)

.PHONY: all pdf figures epub sim-install sim-test hash clean distclean

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

# Provenance manifest: SHA-256 of the book, every figure, and every source/sim file,
# with a UTC timestamp and the current git commit. Commit MANIFEST.sha256 to anchor
# the exact bytes to a point in time (and optionally to a blockchain/notary service).
hash:
	@echo "# Dawn of Eternity -- provenance manifest" > MANIFEST.sha256
	@echo "# generated: $$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> MANIFEST.sha256
	@echo "# git commit: $$(git rev-parse HEAD 2>/dev/null || echo none)" >> MANIFEST.sha256
	@echo "# author: Michael Gronager" >> MANIFEST.sha256
	@find book.pdf book.tex CLAUDE.md README.md CITATION.cff LICENSE.md \
		frontmatter chapters appendices style bibliography \
		figures/scripts figures/pdf sims/src sims/tests notes \
		-type f \( -name '*.tex' -o -name '*.py' -o -name '*.pdf' -o -name '*.bib' \
		-o -name '*.sty' -o -name '*.md' -o -name '*.cff' \) 2>/dev/null \
		| LC_ALL=C sort | xargs sha256sum >> MANIFEST.sha256
	@echo "wrote MANIFEST.sha256 ($$(grep -c '^[0-9a-f]' MANIFEST.sha256) files hashed)"
	@sha256sum MANIFEST.sha256 | sed 's/^/  root hash: /'

clean:
	latexmk -c

distclean:
	latexmk -C
	rm -rf $(BUILD) $(VENV)
