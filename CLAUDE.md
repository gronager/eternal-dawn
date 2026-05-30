# Supraverse — Project Orientation

## What this is

A theoretical cosmology research program built around a single principle: **physics is continuous and conservative, so singularities and discontinuities must be artifacts of incomplete theory, not real features of nature.** Applied carefully, this principle plus Einstein-Cartan gravity plus an infinite quantum vacuum produces a complete cosmological framework with no fine-tuning.

The framework is referred to as **Supraverse Cartasis Theory** (SCT). Cartasis is the bounce membrane — the surface in spacetime where torsion-mediated repulsion overwhelms gravitational attraction and the geometry transitions between parent and child universes.

## What we’re trying to do

Three things, in order of decreasing tractability:

1. **Write the framework down clearly.** Derive each conclusion from minimal axioms. Identify what’s established physics, what’s natural extension, what’s speculative. Get to a state where the framework can be evaluated by people who haven’t been in the conversation that generated it.
1. **Identify computational and observational tests.** Each piece of the framework should make at least one prediction distinguishable from ΛCDM in principle. Catalog these.
1. **Build simulation infrastructure.** Start with 1D spherical Einstein-Cartan collapse and bounce. Progress to perturbation evolution through bounces. Eventually toward supraverse population statistics. The visualization piece (gravity-scaled conformal mapping of the foam manifold) is a parallel track.

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
- LaTeX for the final write-up, markdown for working drafts.

## Repo structure

```
supraverse/
├── CLAUDE.md                    # this file
├── README.md                    # public-facing description
├── 00-axioms.md                 # minimal axioms and what follows from them
├── 01-cartasis.md               # the bounce membrane and its physics
├── 02-baby-universes.md         # how baby universes form and evolve
├── 03-chirality-flip-filter.md  # matter-antimatter dynamics across bounces
├── 04-supraverse.md             # the global structure: foam, equilibrium, lineages
├── 05-dark-sector.md            # dark matter, dark energy, CMB-as-Hawking
├── 06-time-arrow.md             # emergence of time direction
├── 07-observational-tests.md    # predictions and how to check them
├── 08-simulation-plan.md        # tiered simulation pathway
├── 09-visualization.md          # gravity-scaled conformal mapping
├── 10-open-questions.md         # what we don't know
├── notes/                       # working notes, scratch derivations
└── sims/                        # simulation code (eventually)
```

Each numbered file is self-contained enough to be read independently but assumes the axioms in 00. Cross-references use the numbered prefix.

## Workflow notes for Claude Code sessions

When picking up this project:

1. Read CLAUDE.md (this file) and 00-axioms.md first to load context.
1. The conversation that originated this framework was long and rambling. The MD files are the distilled version. Trust the MD files over any prior conversation memory.
1. When extending or revising, preserve the parsimony. The framework’s strength is that it derives a lot from very little. Adding postulates is a regression.
1. When uncertain whether something is established physics vs. speculation, default to labeling speculation. Honesty about epistemic status is more valuable than apparent completeness.
1. LaTeX conversion comes after the markdown is firm. Don’t switch formats prematurely.
1. Simulation code should be Python or Julia, with results saved as plots in `sims/output/`. Numerical relativity libraries (Einstein Toolkit, GRChombo) are too heavy for early exploration; start with custom finite-difference codes.

## What this project is not

- Not a refutation of ΛCDM. ΛCDM works empirically for what it describes. SCT is an alternative that explains a wider set of phenomena from fewer postulates.
- Not a claim that the framework is correct. It’s a coherent worldview that could be right and is testable. The work is finding out.
- Not a vehicle for general philosophical speculation. Each claim should be derivable, computable, or observable in principle.
- Not for popular audiences in its current form. The framework needs to be defensible to professional cosmologists first; popularization comes later if at all.