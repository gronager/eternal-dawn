# Supraverse Cartasis Theory

A working-draft cosmological framework built from minimal axioms.

## The one-paragraph version

Take quantum mechanics, special relativity, and Einstein-Cartan gravity (the natural extension of general relativity that allows torsion to couple to fermion spin). Add an infinite quantum vacuum. The vacuum is unstable: rare large fluctuations reach a critical density (Cartasis density, ~10⁵⁰ kg/m³) where torsion-mediated repulsion overwhelms gravity, and the geometry bounces rather than producing a singularity. The bounce creates a baby universe connected to its parent region through a 3D hypersurface (the Cartasis membrane). The supraverse — the total 4D manifold — condenses into a foam of universes through this mechanism, with equal statistical numbers of matter-biased and antimatter-biased lineages (CPT symmetry preserved globally). Our universe is one bubble in this foam. The CMB is our parent black hole’s Hawking radiation, redshifted by our internal expansion. Dark matter is the gravitational influence of our parent’s matter content (and rejected antimatter content) on us through the Cartasis membrane. The arrow of time is the thermodynamic gradient from our bounce (low entropy) toward our eventual Hawking-mediated dissolution back into the parent’s exterior.

## Why this might be worth taking seriously

The framework derives, from minimal axioms, explanations for:

- The Big Bang (a bounce, viewed from inside)
- Inflation (rapid torsion-driven expansion at the bounce)
- Matter-antimatter asymmetry (statistical inheritance from the seeding vacuum fluctuation, propagated by bounce dynamics)
- Dark matter (parent’s gravitational signature through the membrane)
- Dark energy (parent’s accretion history bleeding through)
- The CMB (parent’s Hawking radiation)
- The arrow of time (low-entropy boundary at bounce, high-entropy future)
- Apparent fine-tuning (typical parameters for observer-supporting universes in the supraverse distribution)
- Information preservation (entanglement across bounce, unitarity maintained globally)

No singularity, no fine-tuned initial conditions, no exotic new fields, no special starting point. Just continuous physics applied carefully.

## Status

This is a research program, not an established theory. It synthesizes work by Poplawski, Penrose, Smolin, Wiltshire, and others, plus original observations about CMB-Hawking identification and supraverse thermodynamics. The framework makes specific testable predictions. None have been rigorously checked yet. The work is in identifying tests, building simulations, and comparing to observational data.

## Where to read more

The monograph is written in LaTeX. Start with `chapters/01-axioms.tex`; each
subsequent chapter develops one piece of the framework.
`chapters/08-observational-tests.tex` lists what would distinguish this framework
from ΛCDM, and `chapters/09-simulation-plan.tex` lays out a tiered computational
pathway. `appendices/A-open-questions.tex` catalogs what we don't yet know.

The pre-LaTeX markdown drafts are kept in `drafts/` for provenance.

## Building

Requires a TeX distribution with `latexmk` (e.g. `texlive-latex-extra`,
`texlive-science`, `lmodern`).

```bash
make          # build book.pdf via latexmk
make clean    # remove LaTeX aux files
make distclean # also remove build/ and the sims venv
```

Layout: `book.tex` is the master document; `frontmatter/`, `chapters/`, and
`appendices/` hold the content; `style/cartasis.sty` is the preamble;
`bibliography/references.bib` holds citations; figures live in `figures/`.