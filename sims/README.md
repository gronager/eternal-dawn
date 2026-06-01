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

# optional: the precision CMB pass (Boltzmann), ~2 s/spectrum, no special hardware
uv pip install -e ".[camb]"               # or: uv pip install camb
.venv/bin/python ../figures/scripts/ch09_camb.py
```

The `camb` extra is optional; the CAMB tests skip cleanly when it is absent, so the
core suite stays dependency-light.

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
   `ρ_C` instead of reaching a singularity (`bounce.py`). Then spherical
   (Weyssenhoff-fluid) collapse matched to an exterior Schwarzschild geometry
   (`os_collapse.py`): a black hole from outside, an expanding cosmology inside.
   - **Tier 1c — the extruder (`extruder.py`).** The bounce *entropy budget*:
     `η_baby = η_parent / D` with `D = S_after/S_before`. A bulk-viscous Eckart
     estimate gives `ln D = 9 ζ̃ (Ω_bounce/T_bounce) J`, `J ≈ 1.38`. The
     prefactor `Ω/T ≈ 1.5e-11` (the same small parameter as the CVE) makes the
     bounce adiabatic to ~1 part in 1e11 for any physical viscosity → `D ≈ 1` →
     descendants **inherit** their parent's ~1 ppb asymmetry; horizon-saturating
     `D ~ 1e49` needs an absurd `ζ̃ ~ 1e12`. *Universes are not born, they are
     extruded.*
2. **Tier 2 — perturbations through the bounce** → primordial spectrum / CMB.
3. **Tier 3 — chirality** → η (CVE + inheritance/dilution, `cve_filter.py`) and
   the dark-matter ratio.
4. **Tier 4 — supraverse population statistics.**

## Reproducibility

Version-controlled, parameter files alongside outputs, validated against analytic
limits, convergence-checked. Outputs in `output/` (gitignored); commit figure
PDFs under `../figures/pdf/` so the book builds without re-running Python.
