# Lattice strong-sector pipeline (Eternal Dawn, Part III / Appendix B)

This is the GPU half of the particle-sector programme. The **analysis** half (turning
measurements into σ, T_c, γ_m, w0) lives in `sims/src/cartasis_sims/lattice.py` and is
unit-tested on CPU (`sims/tests/test_lattice.py`) with synthetic data, so the extraction logic
is validated **independent of any GPU run**. This directory holds the build/run scripts that
produce the configurations and measurements on the GH200.

**Division of labor.** The repo authors and tests the pipeline; you run it on the GH200 and
commit the outputs/logs back; analysis + iteration happen through the repo. Nothing here needs
to be re-derived — it is `make`-and-run, the same ethos as the rest of the book.

---

## The target, and the run order

The gate measurement is **γ_m** (the mass anomalous dimension): it decides whether the candidate
strong sector *walks*, which is the prerequisite for both a small electroweak `S` (Part III) and
the steep mass hierarchy (Part II). But the *run order* starts cheaper, to validate the pipeline
and bank a result before the expensive dynamical-fermion runs:

| # | target | what runs | cost (1× GH200) | module function |
|---|---|---|---|---|
| 1 | **σ** (confinement / string tension) | pure-gauge HMC + Wilson-loop static potential | hours | `static_potential_cornell` |
| 2 | **T_c** (deconfinement = "the condensate melts") | finite-T pure gauge, Polyakov scan | hours–day | `deconfinement_beta_c` |
| 3 | **w0** (scale setting) | gradient flow on the above | minutes | `gradient_flow_w0` |
| 4 | **γ_m** (the GATE) | dynamical-fermion HMC + Dirac mode number | days–weeks | `anomalous_dimension_from_mode_number` |

Steps 1–3 are pure-gauge, single-GPU, interconnect-free — ideal for validating the stack and the
analysis. Step 4 is the prize and the only expensive one.

## The candidate theory

```
gauge group : SU(3)
fermions    : N_f = 2 Dirac in the sextet (2-index symmetric) representation
```
**Why this one:** it is the leading *near-conformal / walking* candidate with a light
flavour-singlet scalar — exactly the composite-Higgs/dilaton Part III needs; SU(3) matches the
colour the Pauli label forces (Chapter 8); and there is substantial existing lattice literature
(Fodor–Holland–Kuti–Nógrádi–Schröder and others) to validate the stack against before trusting an
ED-specific number. **Alternatives** to try if this one disappoints: SU(3) with N_f = 8
fundamental (LatKMI), or minimal walking technicolor SU(2) with N_f = 2 adjoint. The choice of the
walking sector's group/representation is itself part of what is owed (Appendix B); fixing it is a
publishable contribution.

**The gate, concretely:** measure γ_m from the Dirac mode number ν(M) ∝ M^{4/(1+γ_m)}
(Giusti–Lüscher; the stochastic estimator, so no eigenvectors are stored — memory-light,
single-GPU). A clearly nonzero γ_m (order ~0.3, debated in the literature) means the sector walks
and the framework is favourably placed; a vanishing γ_m means it sits in the conformal window or
runs like QCD, and Part III retreats to Part II.

## Hardware notes (GH200)

One GH200 = one Grace-Hopper: 96 GB HBM3 at ~4 TB/s = a strong **single-GPU** node, matched to all
four targets (each lattice fits in one GPU's memory, so **no interconnect is needed**). The build
targets aarch64 (Grace) + Hopper `sm_90`. Lattice solvers are memory-bandwidth-bound and
mixed-precision, so the GH200's bandwidth is the relevant figure of merit, not FP64.

## Build

```bash
lattice/build/build_grid_gh200.sh          # installs deps, builds Grid + its tests, verifies
export GRID=$HOME/ed-lattice/src/Grid/build # the run scripts read $GRID
```

The build script is self-contained: it `apt-get`-installs the external dependencies (GMP, MPFR,
FFTW, HDF5, OpenSSL, autotools), builds **c-lime** from source, clones and configures Grid for
single-GPU Hopper (`sm_90`, no MPI), and crucially runs **`make tests`** — Grid builds only the
library and top-level tests by default, so the sub-directory executables the run scripts call
(`tests/hmc/...`, `tests/core/...`) exist *only* after `make tests`. Step 5 then **lists the test
executables that were actually built**: Grid's test names vary by version, so if a run script
names one that isn't in that list, swap in the matching name. The run scripts call `require_exe`
up front and fail with that same guidance rather than a bare "No such file or directory."

> Grid's stock HMC tests sometimes fix action parameters (e.g. `beta`) in the source rather than
> on the command line. Each run script prints a reminder to check the binary's `--help`; if a
> flag is ignored, edit the `.cc` and rebuild, or supply Grid's HMC XML input.

## Layout

```
lattice/
├── README.md                      # this file
├── build/build_grid_gh200.sh      # deps + Grid + `make tests` + verification (GH200)
└── run/
    ├── _lib.sh                    # shared: $GRID default, require_exe, note_params
    ├── 01_puregauge_potential.sh  # target 1: string tension
    ├── 02_finiteT_polyakov.sh     # target 2: T_c (the melting)
    ├── 03_gradient_flow.sh        # target 3: w0 scale setting
    └── 04_dynamical_modenumber.sh # target 4: gamma_m (the gate)
```

Each run script writes a plain-text/HDF5 measurement file; feed it to the matching
`cartasis_sims.lattice` extractor (or the `--analyze` stub each script prints) to get the number.
