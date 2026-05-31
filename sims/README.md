# Simulations

Calculation and simulation code for Supraverse Cartasis Theory. Python package
`cartasis_sims` (src layout) plus figure scripts under `../figures/scripts/`.

## Layout

```
sims/
├── pyproject.toml
├── src/cartasis_sims/
│   ├── constants.py      # SI constants + Planck 2018 cosmology
│   └── blackhole.py      # Phase-0 "universe = black-hole interior" identities
├── tests/                # pytest
└── output/               # generated results (gitignored)
```

## Setup and run

```bash
cd sims
uv venv && uv pip install -e ".[dev]"
.venv/bin/pytest -q                       # tests/test_identities.py
.venv/bin/python ../figures/scripts/ch02_bh_first_light.py
```

Or from the project root: `make sim-test` and `make figures`.

## Phase 0 — establishing universe == black-hole interior (see Ch. 9)

`blackhole.py` implements the cheap, decisive consistency checks that anchor the
core thesis, in the Copernican spirit (derive from the world's measured numbers):

| corner | result (Planck 2018) | meaning |
|--------|----------------------|---------|
| A | `R_s/R_H = Ω = 1.000` | flat universe sits at its Schwarzschild radius |
| B | `a*_max = 2/Ω ≈ 2`    | order-unity max spin → realistic spin is sub-extremal Kerr |
| C | `T_H ≈ 1.3e-30 K` (= ½ T_GH) | 30 orders below CMB → CMB ≠ present-horizon radiation |
| D | `f = ω_b/(ω_b+ω_c) ≈ 0.157` | membrane filter pass-fraction if DM = rejected matter |

Running the script prints the table, writes `output/first_light.{txt,json}`, and
renders `../figures/pdf/bh_first_light.pdf`.

## Planned tiers (Ch. 9)

1. **Tier 1 — single-bounce Einstein–Cartan dynamics.** Minimal model: the
   modified Friedmann equation `H² = (8πG/3) ρ (1 − ρ/ρ_C)`, which bounces at
   `ρ_C` instead of reaching a singularity. Then spherical (Weyssenhoff-fluid)
   collapse matched to an exterior Schwarzschild geometry: a black hole from
   outside, an expanding cosmology from inside.
2. **Tier 2 — perturbations through the bounce** → primordial spectrum / CMB.
3. **Tier 3 — chirality (flip/filter)** → η and dark-matter ratio.
4. **Tier 4 — supraverse population statistics.**

## Reproducibility

Version-controlled, parameter files alongside outputs, validated against analytic
limits, convergence-checked. Outputs in `output/` (gitignored); commit figure
PDFs under `../figures/pdf/` so the book builds without re-running Python.
