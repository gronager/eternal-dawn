# Simulations

This directory holds simulation code and outputs.

## Structure

```
sims/
├── tier1-single-bounce/        # 1D spherical Einstein-Cartan collapse
├── tier2-perturbations/        # Perturbation evolution through bounce
├── tier3-chirality/            # Flip + filter dynamics
├── tier4-supraverse/           # Population statistics
├── visualization/              # Conformal mapping and rendering
└── output/                     # Generated plots, data files
```

## Getting started

For Tier 1 (first project), start with:

```bash
cd tier1-single-bounce
python -m venv venv
source venv/bin/activate
pip install numpy scipy matplotlib
```

Then implement the spherical bounce as described in `../08-simulation-plan.md`.

## References

Key papers to read before implementing:

- Poplawski, “Cosmology with torsion” series
- Magueijo et al., “Bouncing universes” papers
- Bojowald, loop quantum cosmology bounce calculations (for perturbation analysis)
- Einstein Toolkit documentation (for general numerical relativity context)

## Reproducibility requirements

All simulations should:

1. Be in version control from the start
1. Save inputs (parameter files) alongside outputs
1. Be tested against analytic limits where applicable
1. Pass numerical convergence checks (refine grid, verify stability)
1. Generate plots with matplotlib/Makie that are saved as both raw data (NPZ/HDF5) and figures (PNG/PDF)