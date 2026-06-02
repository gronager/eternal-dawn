# Eternal Dawn — Project Orientation

## What this is

A theoretical cosmology research program built around a single principle: **physics is continuous and conservative, so singularities and discontinuities must be artifacts of incomplete theory, not real features of nature.** Applied carefully, this principle plus Einstein-Cartan gravity plus an infinite quantum vacuum produces a complete cosmological framework with no fine-tuning.

The book is titled **Eternal Dawn**. (The older working title "Supraverse Cartasis Theory" is retired; the `\SCT` macro is being removed.) The bounce membrane — the surface where torsion-mediated repulsion overwhelms gravity and the geometry transitions between a parent region and a child universe — was historically called the *Cartasis* membrane, a pun on **Cartan + κάθαρσις** (catharsis: the bounce dissolves and purifies everything, then extrudes it clean). That pun is now spent **once**, in the bounce chapter; elsewhere we say "bounce density" and "bounce membrane." See **Naming policy** below.

## Central thesis (carry this as the tone, everywhere)

**The dawn is ongoing, not an event.** This is the single largest break from the standard model, and it must color the whole book — not appear as the occasional paragraph. There is no one-time Big Bang finished in the first $10^{-21}$ second. Our bounce membrane is a *persistent* surface at our past horizon, and the parent is *still* extruding matter across it — our "Big Bang" is happening *now*, and goes on until the parent stops feeding (and eventually evaporates). **The dark sector is the live, measurable proof of this ongoing dawn:** dark matter is the parent's mass still pressing through the membrane, dark energy the ongoing *rate* of that inflow. Only an old OGU — parentless, its feed exhausted — has a finished beginning, and (lacking a parent) it has no dark sector. So **dark matter + dark energy ARE the proof of the eternal dawn**, and ΛCDM's finished Big Bang is precisely what we diverge from. Whenever "the Big Bang," "the beginning," or "our dawn" is mentioned, write it as ongoing; reserve the finished, one-shot beginning for the OGU case. The book's title, *Eternal Dawn*, names this.

## What we’re trying to do

Three things, in order of decreasing tractability:

1. **Write the framework down clearly.** Derive each conclusion from minimal axioms. Identify what’s established physics, what’s natural extension, what’s speculative. Get to a state where the framework can be evaluated by people who haven’t been in the conversation that generated it.
1. **Identify computational and observational tests.** Each piece of the framework should make at least one prediction distinguishable from ΛCDM in principle. Catalog these.
1. **Build simulation infrastructure.** Start with 1D spherical Einstein-Cartan collapse and bounce. Progress to perturbation evolution through bounces. Eventually toward supraverse population statistics. The visualization piece (gravity-scaled conformal mapping of the supraverse manifold) is a parallel track.

## Authorship and scope

This is a personal research program, written for my own clarity first and possibly for eventual publication. The intellectual debts include:

- Nikodem Poplawski (torsion bounces, baby universes)
- Roger Penrose (gravitational entropy, conformal cyclic cosmology, Diósi-Penrose collapse)
- Lee Smolin (cosmological natural selection)
- David Wiltshire (inhomogeneous cosmology, timescape model)
- Lior Shamir (JWST galaxy rotation asymmetry observations)
- The standard QFT + GR + Einstein-Cartan literature

The framework is a synthesis, not pure invention. Original contributions are (a) the systematic application of the “no singularities, no discontinuities” principle across measurement, bounces, baryogenesis, dark sector, and time arrow; (b) the identification of CMB with parent Hawking radiation as a direct observable consequence; (c) the gravity-scaled conformal mapping approach to supraverse visualization; (d) the explicit treatment of the supraverse as a thermodynamic equilibrium with universes as entropy-maximizing condensations.

## Style notes

- Direct prose, no hedging where not earned.
- First-principles derivation over numerical fits.
- When something is speculative, label it. When something is established, cite it.
- Plots and simulations should be reproducible from the repo.
- The monograph is authored in LaTeX. Equations and formulas belong in LaTeX, not markdown — author new material directly in the `.tex` chapters. Quick scratch notes may start as markdown in `notes/`, but anything with math should move into the LaTeX chapters promptly.

## Restructure in progress (read this before editing chapters)

The monograph is mid-**reorganization** from the old "as-discovered" order into a clean linear arc that mirrors the plain-language overture. The new manuscript is **`doe.tex` + `doe/`** (canonical-in-progress); the old **`book.tex` + `chapters/`** still compiles and is the source being migrated out of, then deleted. Both build side by side (`make doe` and `make`). Migrate one chapter at a time, building `doe.pdf` green and keeping the 225 sim tests passing. Each `doe/` chapter file carries its **migration contract** as header comments (which old sections feed it, the dedup rule). **Edit `doe/` going forward; treat `chapters/` as the read-only source.**

New structure (see `doe.tex`):

```
0.  First When There's Nothing            (overture — plain language)
PART I  — The Framework
  1. Axioms
  2. No Singularities in Nature           (the bounce; the Cartasis=katharsis pun)
  3. The Void                             (eigenstate, instanton, Boltzmann/Penrose)
  4. Eternal Dawn                         (first universes; spin, baryon, arrow)
  5. Home                                 (child universe; CMB, dark matter/energy)
PART II — Observations and Simulations
  6. What We Can Already See              (live results)
  7. What the Next Decade Will Decide     (near-future)
  8. The Simulation Program
  9. Seeing almost Nothing                (visualization; the void is mostly 3D, tiny 4D pockets)
 10. The Scoreboard                       (vs ΛCDM / Standard Model)
```

## Repo structure

```
cartasis/
├── CLAUDE.md                       # this file
├── README.md                       # public-facing description (opens with the plain-language story)
├── Makefile                        # `make doe` builds the new book; `make` builds the legacy book
├── doe.tex                         # NEW main document (Eternal Dawn) — CANONICAL, in progress
├── doe/                            # NEW chapters (00-overture … 10-scoreboard), Part I + Part II
├── book.tex, chapters/             # LEGACY manuscript — being migrated into doe/, deleted when done
├── frontmatter/                    # title, dedication, epigraph, preface, cover (.tex)
├── appendices/A-open-questions.tex # what we don't know (re-wired into doe after migration)
├── style/cartasis.sty              # preamble: packages, axiom/conjecture envs, macros
├── bibliography/references.bib     # all citations
├── figures/{scripts,pdf,data}/     # figure scripts + generated PDFs/PNGs (shared; do not move)
├── data/                           # observational datasets
├── drafts/, notes/                 # markdown drafts (archival) and scratch notes
├── sims/                           # simulation code (src + tests; reproducible, 225 passing)
└── build/                          # LaTeX build output (gitignored)
```

The new chapters read linearly (each leans on the prior). Cross-references use LaTeX labels (`chap:*`, `ax:*`, `sec:*`, `app:openq`) via `\ref`/`\Cref`. The `drafts/` markdown files are the pre-LaTeX originals — archival only.

## Naming policy

- The book is **Eternal Dawn**; do not reintroduce "Supraverse Cartasis Theory" or `\SCT`.
- Lead with the **warm names** — "first universe(s)," "child universe," "the bounce," "the void." Introduce the formal acronyms **OGU** (original-generation universe) and **BHU** (black-hole universe) once each, then use sparingly.
- "**Cartasis**" appears once, as the Cartan + *katharsis* pun in Chapter 2. Elsewhere: "**bounce density**" (the symbol `\rhoc`, ρ_C) and "**bounce membrane**."

### Migration gardening TODO (do per chapter, and a final sweep)

When migrating each chapter into `doe/`, **garden as you go** — do not just move text:
- **Naming:** purge stray "Cartasis" (keep only the Ch2 pun), `\SCT`, "Supraverse Cartasis Theory"; convert to "bounce density/membrane," "Eternal Dawn," warm names. (Copied stubs like the axioms chapter inherited the old naming — `doe/01-axioms` is gardened; recheck others.)
- **No lab-log voice.** This is a linear monograph, not a research diary. Cut "earlier drafts assumed…," "it is tempting to…," "the honest result is…," "how the sausage is made" framing. State the result and its derivation directly; do not narrate the path that reached it.
- **Figures:** Part I was figure-light early (overture/axioms are prose; the bounce chapter now carries `fig:bounce`; Void has `fig:instanton`; Eternal Dawn has `fig:ogudist`, `fig:cycle`). Part II (observations/sims/viz) is figure-dense. Place each figure in the chapter that *introduces* its idea; reference (don't duplicate) it elsewhere.
- Run a final `grep -rn "Cartasis\|\\SCT\|earlier draft\|tempting" doe/` sweep before deleting the legacy tree.

## Workflow notes for Claude Code sessions

When picking up this project:

1. Read CLAUDE.md (this file), then `doe/00-overture.tex` (the whole story in plain words) and `doe/01-axioms.tex`, to load context. Mind the **Restructure in progress** note above.
1. The conversation that originated this framework was long and rambling. The chapters are the distilled version. Trust the chapters over any prior conversation memory.
1. When extending or revising, preserve the parsimony. The framework’s strength is that it derives a lot from very little. Adding postulates is a regression.
1. When uncertain whether something is established physics vs. speculation, default to labeling speculation. Honesty about epistemic status is more valuable than apparent completeness.
1. The monograph is in LaTeX. Build the new book with `make doe` (`latexmk -pdf doe.tex`); the legacy book still builds with `make`. Author new material in `doe/*.tex`; `chapters/*.tex` is the read-only migration source and `drafts/` is archival.
1. Simulation code should be Python or Julia, with results saved as plots in `sims/output/`. Numerical relativity libraries (Einstein Toolkit, GRChombo) are too heavy for early exploration; start with custom finite-difference codes.

## What this project is not

- Not a refutation of ΛCDM. ΛCDM works empirically for what it describes. Eternal Dawn is an alternative that explains a wider set of phenomena from fewer postulates.
- Not a claim that the framework is correct. It’s a coherent worldview that could be right and is testable. The work is finding out.
- Not a vehicle for general philosophical speculation. Each claim should be derivable, computable, or observable in principle.
- Not for popular audiences in its current form. The framework needs to be defensible to professional cosmologists first; popularization comes later if at all.