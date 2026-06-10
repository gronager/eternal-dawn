# Lattice strong-sector pipeline (Eternal Dawn, Part III / Appendix B)

This is the GPU half of the particle-sector programme. The **analysis** half (turning
measurements into σ, r₀, the spectrum, the bag s_T, f_π) lives in
`sims/src/cartasis_sims/lattice.py` and is unit-tested on CPU (`sims/tests/test_lattice.py`) with
synthetic data, so the extraction logic is validated **independent of any GPU run**. This directory
holds the build/run scripts that produce the configurations and measurements on the GH200.

**Division of labor.** The repo authors and tests the pipeline; you run it on the GH200 and
commit the outputs/logs back; analysis + iteration happen through the repo. Nothing here needs
to be re-derived — it is `make`-and-run, the same ethos as the rest of the book.

---

## The object, and the run order

The object is the **torsiton** — the SU(3)-fundamental baryon (colour from the Pauli label,
Chapter "The Forces from the Field"): three quarks, one of each colour, the chiral soliton the
framework identifies as the elementary fermion. The headline questions are the lattice targets
**L1–L5** of Appendix B: does it confine and bind, and how sharp is the mass-giving condensate
**bag** `s_T` (which sets the whole generation hierarchy). The run order starts cheap, to validate
the pipeline and set the scale, before the dynamical-fermion runs:

| # | target | what runs | run | module function |
|---|---|---|---|---|
| L1 | **σ, r₀** (confinement / scale) | pure-gauge HMC + Wilson-loop static potential | `01` | `static_potential_cornell` |
| — | **T_c** (deconfinement = "the condensate melts") | finite-T pure gauge, Polyakov scan | `02` | `deconfinement_beta_c` |
| — | **w0** (scale setting) | gradient flow | `03` | `gradient_flow_w0` |
| L4 | **the spectrum** (m_π, m_N — the torsiton binds) | dynamical Wilson HMC + baryon/pion correlators | `06`,`07` | `baryon_spectrum`, `correlator_mass` |
| L4 | **the bag s_T** (the generation hierarchy) | dressed-quark profile + condensate 3-pt | `09`,`10` | `bag_profile`, `condensate_3pt` |
| L4 | **generations** (the radial tower count) | variational (GEVP) baryon, nodal basis | `08` | `gevp_spectrum` |
| L2 | **σ = 2πv²** (one substrate, mass + confinement) | pion decay constant f_π | `11` | `decay_constant`, `sigma_2piv2_check` |

Runs `01`–`03` are pure-gauge, single-GPU, interconnect-free — ideal for validating the stack.
`07` is the dynamical ensemble everything downstream measures on.

## The theory

```
gauge group : SU(3)
fermions    : N_f = 2 Dirac, FUNDAMENTAL representation (Wilson), dynamical
```
The torsiton is the fundamental-rep baryon — the object QCD calls the nucleon. SU(3) matches the
colour the Pauli label forces (Chapter "The Forces from the Field"). The strong IR sector is this
torsion-bound soliton, **not** a separate walking gauge theory; the framework's one torsion force
wears the strong and the electroweak faces both. The absolute scale Λ (and with it whether the
propagating-torsion coupling walks) is the **L5** β-function — a separate, harder programme
(Appendix B), not part of this campaign.

**The headline question, concretely:** the fermion-mass hierarchy reduces to one number, the
sharpness `s_T = R0/r0` of the dynamically-massive bag the radial rungs couple to. The lever's
productive window — where the consecutive ladder reproduces the observed ×3477 charged-lepton span —
is `s_T ∈ [0.43, 0.70] r0`. The campaign measures `s_T` (runs `09`, `10`) and tests the sharp
prediction `σ = 2πv²` (run `11`). Results: `RESULTS.md`.

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
single-GPU Hopper (`sm_90`, no MPI), and builds **only the five executables the run scripts need**
— each independently. (Grid's `make` builds just the library + top-level tests, and its full
`make tests` is all-or-nothing: one sibling test that won't compile on a GPU config aborts the
whole suite. Building the five we use one at a time sidesteps that.) Step 5 then **lists the test
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
├── build/build_grid_gh200.sh      # deps + Grid + the 5 test exes + our programs (GH200)
├── src/                           # custom Grid measurement programs (Grid ships no turnkey ones)
│   ├── Makefile                   # compiles against libGrid via grid-config
│   ├── measure_potential.cc       # static potential W(R,T) -> sigma, r0 (L1)
│   ├── generate_dynamical.cc      # dynamical Nf=2 fundamental Wilson HMC (the torsiton vacuum)
│   ├── measure_baryon.cc          # pion + baryon (torsiton) correlators -> the spectrum (L4)
│   ├── measure_baryon_gevp.cc     # variational baryon, nodal basis -> generations (L4)
│   ├── measure_bag_profile.cc     # one-body dressed-quark bag rho(r) -> s_T proxy (L4)
│   ├── measure_condensate_3pt.cc  # genuine 3-body condensate <N|qbar q|N> -> the bag s_T (L4)
│   └── measure_decay_constant.cc  # pion decay constant f_pi -> sigma = 2 pi v^2 test (L2)
└── run/
    ├── _lib.sh                    # shared: $GRID default, require_exe, note_params
    ├── 01_puregauge_potential.sh  # L1: string tension sigma, scale r0
    ├── 02_finiteT_polyakov.sh     # T_c (the melting)
    ├── 03_gradient_flow.sh        # w0 scale setting
    ├── 06_baryon_spectrum.sh      # L4: m_pi, m_N (the torsiton binds)
    ├── 07_dynamical_torsiton.sh   # the dynamical Nf=2 fundamental ensemble
    ├── 08_torsiton_gevp.sh        # L4: generations (variational, nodal basis)
    ├── 09_bag_profile.sh          # L4: one-body bag proxy
    ├── 10_condensate_3pt.sh       # L4: the genuine three-body bag s_T
    └── 11_decay_constant.sh       # L2: f_pi, the sigma = 2 pi v^2 test
```

Each run script writes a plain-text/HDF5 measurement file; feed it to the matching
`cartasis_sims.lattice` extractor (or the `--analyze` stub each script prints) to get the number.
