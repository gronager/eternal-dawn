# SCT status: shown, inferred, open, next

A living summary of where Supraverse Cartasis Theory stands computationally. Honest
about epistemic status: most results are **internal consistency** (the framework
computes a definite, parsimonious answer), a few are **observationally anchored**
(touch data), and the rest are **open**. Every claim below points to the code that
produces it (`sims/src/cartasis_sims/`, figures in `figures/scripts/`).

---

## 1. What we have SHOWN (computed and validated in code)

### The bounce and the black-hole-interior thesis
- **Non-singular Einstein–Cartan bounce** at `ρ_C`, robust across the equation of
  state; collapse → bounce → re-expansion, matched to an exterior Schwarzschild
  (Oppenheimer–Snyder with a torsion bounce). `bounce.py`, `os_collapse.py`.
- **Phase-0 identities** (observationally anchored): `R_s/R_H = Ω ≈ 1`, max spin
  `a* ≈ 2`, membrane filter `f ≈ 1/6`, `T_H` 30 orders below the CMB. `blackhole.py`.

### Baryogenesis: inherited, not manufactured
- **The bounce is adiabatic, D = 1.** Both the homogeneous (bulk-viscous) and the
  inhomogeneous (shear-dominated) channels give entropy dilution `D − 1 ≲ 1e−10`,
  because the bounce is far faster than the thermal time (the same `Ω/T ~ 1e−11`
  that limits the CVE). The same torsion that removes the singularity tames the BKL
  entropy divergence. `extruder.py`.
- **η is inherited debris:** `η_child = η_infall = b·η_parent`, the asymmetry of
  what the hole *ate*. Two families — clean/photon-dominated (`η ≪ 1`, like us) vs
  baryon-rich/degenerate (`η → 1`). `accretion.py`.
- **CVE estimate:** spin-sourced asymmetry `η ~ C(ω/T)`, `ω/T ~ 1.5e−11` at maximal
  spin → 41× short *if* the CVE alone had to source η. `cve_filter.py`.

### Population: where we sit
- **We are BHU1 (robust).** Generation-dependent branching: an OGU is super-critical
  (many children) but its children sit at `M_vis` and are sub-critical, so
  `P(BHU_n | n≥1) = (1−ε)ε^{n−1}`, `P(BHU1) ≈ 0.9` for `ε ~ 0.1` — **independent of
  OGU size.** `population.py`.
- **Mass-budget caps branching:** viable children need ~horizon-mass holes, so
  `m = ε·M_parent/M_vis`, not `~1e18`. A narrow viable band ⇒ shallow. `population.py`.

### Growth and OGU genesis
- **Causal rate cap:** the M² runaway self-limits at `dM/dt ≤ c³/2G ~ 3e12 M⊙/yr`,
  so `M ~ c³t/2G`. OGUs (infinite void) grow steadily and never run out; internal
  holes (finite parent) deplete and freeze. `growth.py`.
- **OGU mass = horizon mass = supraverse-age clock:** `M_OGU ~ c³ t/2G`. Our mass at
  ~1 Hubble time; `1e65 kg ↔ ~1e12 Hubble times ↔ ~1e11 siblings`. Vortex spin
  window (survival floor + centrifugal choking) gates genesis. `genesis.py`.
- **Average OGU size = birth rate, via KJMA:** `M_avg ~ (c²/2G)(c/β)^{1/4}`; the void
  tiles into a Johnson–Mehl foam (impinged cells are polygons, not circles).
  OGUs are effectively immortal (`t_H ~ 1e143–1e178 s`). `void_foam.py`.
- **Packed-foam coarsening by radiative siphoning.** Once accretion halts, each OGU
  sits at marginal Hawking balance (bath ~ neighbours' radiation ~ `T_dS ~ 1e-30 K`);
  negative heat capacity makes it unstable, so big (cold) OGUs siphon small (hot)
  ones, which evaporate into them — Ostwald ripening for universe-holes, on the
  Hawking timescale. So OGUs are immortal on any ordinary clock but not eternal.
  `void_foam.py` (`simulate_coarsening`).
- **The recursive cycle (speculative).** A heat-dead OGU's smooth de Sitter far
  future is a low-Weyl-entropy void that nucleates the next foam of OGUs — the void
  is not primordial but the asymptotic state every universe reaches. The apparent
  paradox (heat death = max entropy) is resolved by Penrose: *total* entropy rises
  monotonically, but *Weyl/clumping* entropy cycles low→high→low and resets. SCT is
  "Penrose CCC made to branch." `cycle.py`.
- **An OGU is a mortal recursive timeline.** It grows to the **geometric** cap (the
  Nariai mass `c³/3√3GH_Λ ~ a Hubble mass`, where the BH and cosmological horizons
  merge) in ~1 Hubble time. It is **never in true equilibrium** — Schwarzschild–de
  Sitter has T_b > T_c always (equal only at the unstable Nariai point), so it always
  **evaporates** on `t_evap ~ 1e142 s` (a near-Nariai *lukewarm plateau* mimics
  equilibrium but isn't). This corrected an earlier "thermal equilibrium, growth
  stops" slip. The interior does **not** end when the membrane evaporates: the membrane
  is the baby's *past* horizon, and exterior evaporation is dual to the interior
  reaching its own de Sitter future (an infinite interior time maps to the finite
  exterior horizon lifetime) — the umbilical pinches off, nothing's clock stops. In its
  lifetime an OGU hosts **~1e34 interior cycles**: a *fractal of finite, pulsating
  timelines*, every membrane mortal, the tree eternal, the flat void the one eternal
  substrate. `recursion.py`, `sds.py`.

### The OG birth rate, and the instanton (computed -> dilute supraverse)
- **β = λ·H_Λ⁴/c³.** Via the recursive cycle, a de Sitter nucleation rate; the
  **prefactor is pinned by the observed dark energy** (~4e-97 m⁻³s⁻¹). `birth_rate.py`.
- **The instanton action is computed: I ~ 2.5e84.** OG nucleation is an upward de
  Sitter fluctuation, so `I = 2π M_seed c²/ℏH_Λ` (Boltzmann = horizon-entropy
  deficit), with the minimal self-gravitating bounce seed `M_seed ~ 9e14 kg`. This is
  **84 orders above the percolation threshold I_crit ~ 1.4** (robust). So `λ ≈ 0` and
  **the supraverse is DILUTE** — isolated, round, eternal-inflation-like island
  universes, NOT a packed foam. The packed-foam constructions (avg size, coarsening,
  polygonal render) are the non-realised percolating limit; BHU1, the recursive
  cycle, and `M_OGU ~ c³t/2G` survive. `birth_rate.py` (`seed_instanton_action`).

### CMB (Tier 2) — observationally anchored
- **Scale-invariant primordial spectrum from the bounce.** Mode evolution through
  the bounce gives `n_T = 3 − |2p−1|`, `p = 2/(1+3w)`; matter contraction (w=0)
  ⇒ `n_T = 0` **exactly**, with no inflaton (validated <1% vs analytic). EC supplies
  the non-singular bounce the matter-bounce scenario otherwise lacks. `primordial.py`.
- **Acoustic peak positions reproduced.** Shared recombination physics: `r_s ≈ 142
  Mpc`, peaks at `ℓ ≈ 230, 540, 850, 1150, 1460` vs Planck — a few percent. `acoustic.py`.
- **Acoustic peak heights reproduced (semi-analytic).** Baryon loading `R ≈ 0.62`
  (computed) drives the odd/even alternation; a tight-coupling `D_ℓ` (monopole +
  Doppler + driving + Silk damping) on the bounce's scale-invariant input matches the
  Planck peak heights to tens of percent. The 1st/2nd ratio ≈ 2.2 **measures ω_b** —
  which SCT sources from the inherited η; ω_c-sensitive 3rd/4th peaks need a Boltzmann
  code. `cmb_peaks.py`.

### Dark sector
- Dark energy as parent accretion confronts DESI DR2 (`dark_energy.py`); CMB-axis /
  galaxy-spin alignment as the cheapest decisive test (`axes.py`, `galaxy_spins.py`).
- **Dark matter as look-back: gravity felt before light seen** (`lookback.py`). The
  membrane is our *past* horizon, so parent-injected matter arrives in our deep past.
  The particle horizon grows (47→64 Gly asymptote; +1 Gyr → +1 Gly further back), but
  the parent's mass is already in our metric (felt as dark matter) while its light is
  still climbing in — gravity leads light. The dark sector is a *rate*: dark matter ~
  accumulated projected parent mass, dark energy ~ d/dt of its growth.

### Big picture
- **"From what we observe, we must exist" (epistemic Copernicus).** Ch1 synthesis:
  three measured-on-Earth ingredients (fluctuating vacuum, gravity, spin→torsion) plus
  infinity *force* a recursive supraverse that births universes like ours forever —
  the claim is not "might be" but "must be," and the Copernican principle is returned
  as a result, not assumed.

---

## 2. What we have INFERRED (follows, with stated assumptions)

- **We are BHU1**, one of an OGU's (possibly enormous) brood of siblings.
- **We descend from a fair-sample / horizon-scale progenitor** (clean family), not a
  star-fed one — our low `η` requires `b ~ 1` accretion.
- **η is a lineage *sign* invariant**; its magnitude tracks accretion, not a fixed
  number passed down.
- **The 41× is not a problem** — inheritance threads the asymmetry; the CVE is
  demoted from *source* to *sign-biaser* (it still ties baryon sign to the spin axis).
- **Lineages are chirally pure** (CPT-even bounce; sign set once at the OG seed).
- **The OGU mass, supraverse age, sibling count, and average foam size all reduce to
  one unknown** — the OG birth rate β (equivalently the supraverse age).
- **The foam is immortal and slowly coarsening**; growth decelerates from the
  birth runaway to a merger/expansion crawl, never quite stopping.

---

## 3. KEY OPEN UNKNOWNS (the few numbers everything reduces to)

1. **β — the OG birth rate (Q12).** Sets M_OGU, supraverse age, sibling count, and
   the average foam cell size. The single biggest lever; currently unpinned.
2. **The void density / temperature (Q15c).** Sets the growth time `τ_gen ∝ T⁻³`,
   hence the depth cap and whether the void is condensed (horizon-limited OGUs).
3. **C (anomaly/sphaleron coefficient) and η_min.** Fix the OGU founder's `η` value
   — why `~6e−10` and not `1e−3` or `1e−15`. The genuinely open part of baryogenesis.
4. **ε (viable fraction per universe).** Sets the branching ratio and `P(BHU1)`.
5. **ρ_C precise value and its matter-dependence** (Q1, Q2).

---

## 4. NEXT STEPS

### CMB (highest impact — the decisive test)
- **Boltzmann peak *heights*** (CAMB/CLASS-class): turn "positions ✓" into a full
  C_ℓ curve over Planck data. Sensitive to ω_b, ω_c, and the tilt/amplitude.
- **Scalar spectrum + tensor-to-scalar ratio r** from the bounce (we did tensor);
  pin the slightly-soft (w≲0) contraction that gives n_s = 0.965.

### Mechanism
- **Tier 1 inhomogeneous collapse** carrying its own entropy budget, including
  gravitational particle production at the bounce (the one entropy source not yet in
  the D=1 result).
- **Vortex entrainment from real Kerr–Cartan dynamics** (Q3): does the rotating
  bounce actually pump void matter, and at what efficiency? Would turn M_OGU from
  "open, age-clocked" into a number given the void density.

### Observational
- **The galaxy-spin / CMB-axis alignment test** (Shamir + Planck axes) — the cheapest
  potentially-decisive measurement; settle whether the galaxy-spin axis favours the
  CMB axis over the Galactic pole.

### Population
- **Pin β** from any observable handle on the foam (or invert from a measured OGU
  property), collapsing unknowns 1–2.

---

*Through-line of the recent work: baryogenesis (inherited, adiabatic) → population
(we are BHU1) → growth (causal-rate, horizon-mass) → genesis (vortex + horizon mass)
→ foam (average size = birth rate) → CMB (scale-invariant spectrum + acoustic peaks).
The framework now has an end-to-end, mostly-computed spine; the open parts are a
handful of numbers (β, void density, C) and the CMB peak heights.*
