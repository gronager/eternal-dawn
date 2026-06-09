# The Horizon — particle genesis: design & roadmap

**Goal.** Watch the Standard Model condense out of the cooling chiral fluid: leptons, quarks,
baryons and mesons emerging as *distinct* objects from ONE transition, with relative masses and the
genesis abundance. Classical stochastic effective field theory — the *mechanism and the order of
magnitude*; the lattice (L4) remains the exact arbiter for the quantum spectrum.

## The framework's claim we are realizing

- All fermions are **torsitons** (chiral solitons); **one** condensate, **one** scale Λ.
- `m(tower, gen) = Λ · c_T · |O_n(s_T)|²` — a per-tower coupling `c_T` (colour/charge grip) times a
  per-generation radial rung `O_n(s_T)` (the bag sharpness).
- The strong-torsion sector **confines** (string tension `σ = 2π v²`): coloured quark-knots bind into
  colour-singlet baryons; leptons are colourless and free.
- **One** transition (the unification), unlike the Standard Model's staged electroweak-then-QCD.

## What is BUILT (the foundation — all tested, all runnable)

| module | what it gives | status |
|---|---|---|
| `genesis_quench` (2D), `genesis_quench3d` (3D) | the condensate forming; baryons (Skyrmions) condensing; `M ∝ √κ`; Kibble–Zurek abundance | ✅ |
| `genesis_species` | distinct chiral **sectors** (towers), **batched/GPU**, species masses OUT of the soliton dynamics | ✅ |
| `genesis_confine` | **confinement** — quarks confined, baryons/mesons/leptons free; tied to `σ` | ✅ (3a) |
| `horizon` | the spectrum **factorisation** `c_T × rung`; the twelve-fermion cascade from one transition | ✅ |
| `fermion_masses`, `generations`, `self_consistent`, `chiral_quark_soliton` | the configmass / `s_T` machinery | ✅ |

## What the FULL Step 3 still needs (owed pieces, hardest first)

1. **SU(3)-flavour Skyrme baryon octet/decuplet** — the *real* baryon spectrum (p, n, Λ, Σ, Ξ; Δ, Σ\*,
   Ξ\*, Ω). The chiral field `U ∈ SU(3)`, the B=1 Skyrmion via the SU(2) embedding, then **collective
   quantization** (the SU(3) rotor on the collective coordinates) → the octet/decuplet with mass
   splittings from the strange-quark mass. *Difficulty: HIGH* — the splitting is a quantum rotor
   problem, not a field sim; established (Guadagnini, Witten) but intricate.
2. **SU(3)-colour flux tubes** — generalise the U(1)-colour confinement proxy (3a) to real colour.
   *Difficulty: HIGH* (lattice-grade; the proxy already captures the mechanism).
3. **The electroweak / charge sector** — leptons as colourless knots carrying electric charge (the
   charged-lepton and neutrino towers, `c_T` from the charges). *Difficulty: MEDIUM.*
4. **The couplings from the lattice** — `Λ, σ, s_T, c_T` → absolute masses and the full hierarchy.
   *In progress* (the L4 ensemble).

## Staged execution plan (recommended order)

- **3a ✅ confinement demonstrator** — quarks confined, baryons/leptons free (`genesis_confine`).
- **3b ✅ the baryon octet+decuplet** (`genesis_su3`). The hedgehog soliton + rigid-rotor collective
  quantization gives the nucleon, the Δ, and the hyperons (Λ, Σ, Ξ; Σ*, Ξ*, Ω). N, Δ, Ω calibrate
  (f_π, e, the strange scale); the rest are **predictions** — the decuplet equal spacing to ~1%, the
  octet to ~3–9%. *The Horizon's first genuine prediction.* Known limitation: the Λ–Σ split (their
  isospin) needs the full Yabu–Ando treatment; here they are degenerate. With f_π, e from the lattice
  Λ this becomes absolute-scale from first principles.
- **3c ✅ the lepton sector** (`genesis_leptons`). The colourless half of the soliton: free (3a), no
  flavour rotor, so the spectrum is the RADIAL tower — the generations e, μ, τ from the s_T ladder; the
  neutrinos are the same rungs with the weakest condensate grip → sub-eV, no new scale. The unified
  picture (`horizon_species`): the whole observable fermion spectrum from one soliton — leptons
  (colourless rungs) + baryons (coloured Skyrmions), 11 orders of magnitude, two organising principles
  set by colour.
- **3d ✅ the unified genesis** (`genesis_horizon`). The capstone: all observable matter from ONE
  cooling transition — the colourless leptons (rungs) + neutrinos (3c) and the coloured sector confined
  (3a) into the baryon octet/decuplet (3b), every species switching on together as the condensate
  forms. 14 species, 11 orders of magnitude, one transition. **Parametrised for the lattice:**
  `genesis(scale_MeV=Λ, s_T=s_T)` runs on the calibrated values now and re-runs for the absolute-MeV
  prediction the moment L4 lands — the re-run hook.
- **3e → feed the lattice numbers** (`Λ, σ, s_T, c_T`) → absolute MeV + the genesis abundance. The
  *prediction* — the fan either lands on the observed twelve or it does not.

## The convergence

Everything keys on the L4 ensemble (cooking now): **Λ** (scale → absolute masses), **σ** (string →
confinement), **s_T** (bag → the generation hierarchy), **c_T** (tower couplings). When those land,
3e turns the relative cascade into a falsifiable prediction — and the *same* numbers fix the cosmology
(genesis abundance) and the particle masses at once.

## Honest scope (standing)

- Classical stochastic EFT: formation, masses to ~30%, abundance, topology — **not** the exact quantum
  spectrum (the lattice). The octet (3b) is the established Skyrme result (~30% on the baryon spectrum
  — order of magnitude, the stated bar).
- The full ×3477 fermion hierarchy needs the `s_T` structure; it is too extreme to resolve as a single
  classical knot-size spread in one box — it lives in the configmass/lattice, not the genesis sim.
- Confinement here is a U(1)-colour proxy (exact-in-2D linear string); real SU(3) flux tubes are 3b/3c.

## Compute (the GPU)

The cooling is sequential in *time* (a Langevin Markov chain), but **every per-step operation is
parallel** over space and over the batch. `genesis_species` is already batched and backend-agnostic
(`numpy`↔`cupy`, `CARTASIS_GPU=0` forces CPU). The 96 GB holds huge grids (512³ ≈ 2 GB) and many
batched realizations (species, seeds, Kibble–Zurek rates) at once — the embarrassingly-parallel
dimension that fills the GH200. For production: `pip install cupy-cuda12x`, run, done.

## Next action

When the lattice lands: `run/01` for `r0/a`, then the `run/10` ground-state line for `s_T` and `Λ`.
Then **3b** (the SU(3) octet) is the build to start — it is the first piece that turns the Horizon
from a *picture* into a *prediction*.
