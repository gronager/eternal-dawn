MAIN := book
BUILD := build
SIM := sims
VENV := $(SIM)/.venv
PYTHON := $(VENV)/bin/python
FIGSCRIPTS := $(wildcard figures/scripts/ch*.py)

.PHONY: all pdf doe figures epub sim-install sim-test hash clean distclean

all: pdf

pdf:
	@mkdir -p $(BUILD)
	latexmk -pdf $(MAIN).tex

# Dawn of Eternity restructure (built alongside book.tex during migration).
doe:
	@mkdir -p $(BUILD)
	latexmk -pdf doe.tex

# Regenerate every chapter figure (PDFs land in figures/pdf/).
figures: sim-install
	@for f in $(FIGSCRIPTS); do echo ">> $$f"; $(PYTHON) $$f || exit 1; done

# EPUB3 via pandoc. Pandoc reads the LaTeX directly and embeds the figures as PNGs
# (EPUB cannot carry PDF images), resolving them from figures/pdf/ with the .png
# twins that the figure scripts already produce alongside each .pdf. Install pandoc
# from your package manager (apt/brew 'pandoc') -- a single binary, no TeX needed.
# Run `make figures` first if any figure PNGs are missing.
epub:
	@command -v pandoc >/dev/null 2>&1 || { \
		echo "pandoc not found. Install it (apt/brew 'pandoc') -- a single binary."; \
		exit 1; }
	pandoc $(MAIN).tex --from=latex --to=epub3 \
		--resource-path=.:figures/pdf --default-image-extension=png \
		--mathml --toc --toc-depth=2 --number-sections \
		--top-level-division=chapter \
		--epub-cover-image=figures/pdf/cover.png \
		--metadata title="Dawn of Eternity" \
		--metadata author="Michael Gronager" \
		--metadata lang=en \
		-o $(MAIN).epub
	@echo "wrote $(MAIN).epub ($$(du -h $(MAIN).epub | cut -f1)) -- tracked, downloadable from the repo"

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
	@find book.pdf book.epub book.tex CLAUDE.md README.md CITATION.cff LICENSE.md \
		frontmatter chapters appendices style bibliography \
		figures/scripts figures/pdf sims/src sims/tests notes \
		-type f \( -name '*.tex' -o -name '*.py' -o -name '*.pdf' -o -name '*.epub' \
		-o -name '*.bib' -o -name '*.sty' -o -name '*.md' -o -name '*.cff' \) 2>/dev/null \
		| LC_ALL=C sort | xargs sha256sum >> MANIFEST.sha256
	@echo "wrote MANIFEST.sha256 ($$(grep -c '^[0-9a-f]' MANIFEST.sha256) files hashed)"
	@sha256sum MANIFEST.sha256 | sed 's/^/  root hash: /'

clean:
	latexmk -c

distclean:
	latexmk -C
	rm -rf $(BUILD) $(VENV)
