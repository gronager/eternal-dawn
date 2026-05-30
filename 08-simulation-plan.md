# 08 — Simulation Plan

## Philosophy

The framework’s strength is that it should be simulable from first principles. Rare large vacuum fluctuations + Einstein-Cartan bounces + thermodynamic equilibrium are all computable problems. The path from “interesting framework” to “constrained by data” runs through simulation.

The simulation work is tiered. Each tier is a tractable project of definite scope. Higher tiers build on lower ones. A small group could complete tiers 1-2 in a year; tiers 3-4 are multi-year programs.

## Tier 1: Single-bounce Einstein-Cartan dynamics

### Goal

Numerically evolve a spherical collapse of fermionic matter through the bounce. Output: complete characterization of the bounce as a function of input parameters.

### Setup

- 1D spherical symmetry (radial coordinate $r$, time $t$)
- Initial conditions: uniform ball of matter at moderate density, slightly overdense to drive gravitational collapse
- Matter modeled as Weyssenhoff fluid: perfect fluid with spin density (semi-classical treatment of fermion spin in continuous limit)
- Evolution: Einstein-Cartan field equations with torsion sourced by spin density
- Tracking: density, pressure, spin density, metric components, torsion components as functions of $r$ and $t$

### Numerical method

- Finite difference scheme on a radial grid, ~$10^3$ to $10^4$ points
- Adaptive mesh refinement near the bounce (where gradients become extreme)
- Implicit time stepping for stability near peak density
- Conservation checks: total mass-energy, total spin, total particle number

### Expected output

- The bounce happens at $\rho \sim \rho_C \sim 10^{50}$ kg/m³ (verified)
- Bounce timescale in proper time
- Bounce duration and profile shape
- Maximum density reached as a function of initial conditions
- Equation of state through the bounce
- Energy balance: how much kinetic energy is in expansion vs. contained in spin/torsion fields

### Tools

- Python with NumPy/SciPy for prototyping
- Optionally Julia for performance
- Standard ODE solvers (RK4, implicit methods)
- Custom grid handling and AMR
- Plotting with matplotlib

### Effort estimate

- 2-3 months for a competent grad student
- ~500-1500 lines of code
- Validation against existing single-bounce papers (Brandt, Magueijo, Poplawski)

### Status of existing work

Several authors have done pieces of this. No one has produced a clean, well-documented reference implementation. **A valuable side-output of this tier is producing such a reference implementation.**

## Tier 2: Perturbation evolution through the bounce

### Goal

Take output of Tier 1 (clean bounce profile) and evolve cosmological perturbations through it. Output: primordial power spectrum from the bounce, comparable to CMB observations.

### Setup

- Start with Tier 1 spherical bounce background
- Add linear perturbations to the metric and matter fields
- Evolve perturbations using linearized Einstein-Cartan equations
- Track power spectrum of scalar, vector, tensor modes through the bounce
- Compute predictions for CMB observables: scalar spectral index $n_s$, tensor-to-scalar ratio $r$, running of spectral index, non-Gaussianities

### Numerical method

- Mode-by-mode evolution in Fourier space
- Long-wavelength modes treated as classical fields
- Short-wavelength modes treated quantum-mechanically (mode functions)
- Bunch-Davies-like initial conditions at large negative time (deep in collapse)
- Track modes through bounce and into post-bounce expansion

### Expected output

- Primordial power spectrum amplitude and shape
- Spectral tilt
- Tensor-to-scalar ratio
- Comparison to Planck CMB observations:
  - $n_s = 0.965 \pm 0.004$
  - $r < 0.06$ (upper limit)
- Comparison to ΛCDM + inflation predictions

### Decisive test

If the bounce produces a near-scale-invariant spectrum with $n_s \approx 0.96$ and $r$ within observational bounds, **the framework is observationally viable.** If the spectrum is dramatically different, the framework fails this test.

This is the central calculational challenge for SCT. Existing work on bounce cosmologies (loop quantum cosmology, ekpyrotic models) provides templates for similar calculations.

### Effort estimate

- 6-12 months for someone with cosmological perturbation theory background
- ~1000-3000 lines of code building on Tier 1
- Reference: existing LQC bounce perturbation calculations (Bojowald and others)

## Tier 3: Chirality dynamics through the bounce

### Goal

Quantify the flip + filter mechanism. Predict the chirality bias of a baby universe given the parent’s properties.

### Setup

- Use Tier 1 bounce background
- Include fermion content explicitly (not just Weyssenhoff fluid)
- Track baryon number, lepton number, chiral currents through the bounce
- Include electroweak physics at temperatures above EW scale (sphaleron processes)
- Compute filter efficiency for chirality
- Compute net baryon number production / preservation

### Decisive tests

- Does the framework predict $\eta = n_B/n_\gamma \sim 6 \times 10^{-10}$?
- Does the framework predict dark matter / baryon ratio $\sim 5$?
- Are these predictions robust to parameter choices, or do they require fine-tuning?

### Effort estimate

- 12-24 months
- Significant theoretical work in addition to numerical
- Coupling Einstein-Cartan to Standard Model is non-trivial

## Tier 4: Supraverse population statistics

### Goal

Simulate the full supraverse: many vacuum fluctuations producing OG bounces, each seeding lineages of descendant universes. Compute equilibrium distributions of universe parameters. Check that observer-supporting universes cluster near our observed parameters.

### Setup

- Treat the supraverse as a stochastic process: vacuum fluctuations produce OG bounces at random locations
- Each OG bounce seeds a lineage
- Each universe evolves cosmologically and produces descendants through internal black hole formation
- Track populations of universes across many generations
- Compute distributions of: universe mass, age, baryon asymmetry, chirality lineage depth

### Decisive tests

- Equilibrium distribution exists and is computable
- Observer-supporting universes (those with viable parameters) cluster around some peak
- Our observed parameters sit near the peak (validating fine-tuning resolution)

### Effort estimate

- Multi-year project, several people
- Requires Tier 1-3 to be working as input
- Output: most comprehensive cosmological prediction the framework can make

## Tier 5: Direct observation comparisons

### Goal

Use the simulations to make specific predictions for observational programs (Planck, DESI, Euclid, LiteBIRD, JWST, etc.) and compare directly.

### Examples

- Predict precise CMB spectrum shape, including specific deviations from ΛCDM that should be present
- Predict large-scale structure clustering with dark matter as parent influence
- Predict Hubble parameter at different redshifts
- Predict galaxy rotation asymmetry magnitude from parent rotation
- Predict dark energy evolution trajectory

### Decisive tests

The framework is genuinely tested against existing and forthcoming observational data. Each prediction either matches or doesn’t.

### Effort estimate

Each comparison is its own project. Multi-decade research program at full scope.

## Minimum viable project

If we want one project to do that demonstrates the framework’s viability:

**Tier 2 (perturbation through bounce) with focus on CMB power spectrum.**

This is decisive: either the bounce reproduces something close to the observed CMB spectrum or it doesn’t. The result is publishable in either case. The effort is bounded at ~12 months. The required expertise (numerical relativity + cosmological perturbation theory) is widely available.

This should be the first priority for any simulation effort.

## Visualization parallel track

Alongside the physics simulations, the visualization work (see `09-visualization.md`) can proceed independently. It doesn’t require Tier 2 outputs to start — it can begin with simpler models of supraverse structure and refine as physics simulations produce better inputs.

The visualization track is more accessible and might be the right place for someone with computational/graphics interests to start. It produces “supraverse art” that’s both scientifically meaningful and aesthetically striking.

## Tooling

### Languages

- **Python:** prototyping, plotting, analysis. Default for Tier 1.
- **Julia:** performance-sensitive numerics. Optional for Tier 1-2, advisable for Tier 3-4.
- **C++ with Einstein Toolkit:** if existing numerical relativity infrastructure is leveraged. Heavy but well-tested.

### Libraries

- NumPy, SciPy, matplotlib (Python)
- DifferentialEquations.jl (Julia)
- Einstein Toolkit, GRChombo (numerical relativity)
- CAMB or CLASS (CMB power spectrum, for comparison to standard predictions)

### Hardware

- Tier 1-2: workstation or single-node cluster, ~hours-days runtime
- Tier 3: small cluster, ~days-weeks runtime
- Tier 4: significant cluster resources, ~weeks-months runtime
- Tier 5: depends on specific calculations

### Storage

- Simulations should produce reproducible outputs in standardized formats (HDF5 preferred)
- Visualization data should be portable across plotting backends
- Code in git, with versioning of major releases

## Reproducibility

All simulation code should be:

- Open source from the start
- Documented well enough that others can reproduce results
- Validated against analytic limits where they exist
- Tested for numerical convergence (refine grid, check stability of results)

The point isn’t to “prove” SCT through simulations — the point is to make the framework testable. Reproducible code lets others check, extend, and challenge results.

## Cross-references

- `00-axioms.md`: physical foundation for what’s being simulated
- `01-cartasis.md`: bounce physics (Tier 1)
- `02-baby-universes.md`: baby evolution (Tier 2)
- `03-chirality-flip-filter.md`: chirality dynamics (Tier 3)
- `04-supraverse.md`: population structure (Tier 4)
- `07-observational-tests.md`: targets for Tier 5
- `09-visualization.md`: parallel track