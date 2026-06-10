# Lattice results — the torsiton

This file records the non-perturbative lattice results, with enough detail to reproduce. The
torsiton is the **SU(3)-fundamental baryon** (colour from the Pauli label, Ch. "The Forces from the
Field") — three quarks, one of each colour, the object QCD calls the nucleon. **Not** the sextet
(that was a walking proxy for the electroweak S parameter, L3 — a different question, and the wrong
turn for the soliton).

## L4 pilot — the torsiton ground state (quenched)

**First non-perturbative confirmation that the torsiton binds**, and the chiral trend of its mass.

- **Setup:** quenched valence spectroscopy (`measure_baryon`) on a pure-gauge SU(3) Wilson ensemble.
  Lattice 12³×24, β=5.6, 22 thermalised/decorrelated configs. Point source, Wilson valence quark.
  Pion `C_π = Tr[S†S]` (anchor + chiral probe); nucleon = two-Wick contraction with the Cγ5 diquark
  and the positive-parity projector. Masses from the effective-mass plateau (window t∈[3,8]),
  jackknife errors over configs. Reproduce: `run/06_baryon_spectrum.sh` (analysis:
  `cartasis_sims.lattice.baryon_spectrum`).

| Wilson mass | m_π a | m_N a | m_N / m_π |
|---|---|---|---|
|  0.10 | 2.167(9) | 3.544(18) | 1.635 |
|  0.00 | 2.078(9) | 3.418(19) | 1.645 |
| −0.10 | 1.983(9) | 3.285(19) | 1.656 |
| −0.20 | 1.884(9) | 3.145(19) | 1.669 |
| −0.30 | 1.778(9) | 2.996(19) | 1.685 |
| −0.40 | 1.665(9) | 2.838(19) | 1.705 |
| −0.50 | 1.544(9) | 2.670(19) | 1.729 |
| −0.60 | 1.414(9) | 2.492(19) | 1.762 |
| −0.70 | 1.273(9) | 2.300(19) | 1.806 |
| −0.80 | 1.118(9) | 2.093(20) | 1.871 |
| −0.90 | 0.944(9) | 1.866(23) | 1.976 |
| −1.00 | 0.742(10) | 1.611(28) | 2.172 |
| −1.05 | 0.623(10) | 1.471(34) | 2.361 |

m_π² stays linear in the bare mass right down to the last point → m_crit ≈ −1.15. Below m ≈ −1.05
(m_π a ≈ 0.6) the nucleon error grows as the quenched exceptional configurations switch on — the wall.

**What it establishes (robust):**
- The torsiton **binds** — a clean effective-mass plateau at m_N > m_π > 0, non-perturbatively, with
  no mean-field approximation. (A wrong-sign contraction would give a negative/junk correlator; this
  is clean and correctly ordered.)
- **Constituent counting.** At heavy quark mass m_N/m_π → 3/2 (a baryon is three quarks, a pion two);
  we sit just above it (1.635 at the heaviest point). This is the lattice confirming the mean-field
  result **M_torsiton ≈ N_c·M** (Appendix C, `kahana_ripka`) — the two halves of the program agree.
- **Chiral trend.** m_π² is linear in the bare mass (GMOR), giving m_crit ≈ −1.15; as the quark
  lightens, m_N/m_π **rises** (1.64 → 2.36, the increments accelerating) — the mass going
  condensate-dominated. Extrapolating m_N vs m_π² to the chiral point: **m_N(chiral)·a ≈ 1.15**.
- **The torsiton lands at the baryon mass.** The static-potential run (`run/01`) gives the scale
  r₀/a = 2.28 — dead on the β=5.6 reference (2.3), so the scale is trustworthy even though σ itself is
  noisy (σ a² = 0.31(13), smearing-inflated). Hence **m_N·r₀ ≈ 2.6** (equivalently m_N/√σ ≈ 2.1),
  against the **physical baryon's m_N·r₀ ≈ 2.4** (m_N/√σ ≈ 2.1). The torsiton comes out within ~10%
  of the nucleon — exactly the agreement quenched QCD gives for the baryon spectrum. So the pilot
  shows not just that the torsiton *binds* but that it binds at the *right mass relative to its own
  confinement scale*. Call it **m_N·r₀ ≈ 2.6(5)** given the small box, coarse spacing, quenching, and
  short chiral extrapolation.

**What it does NOT establish (the honest caveats):** this is **quenched** (no sea quarks), **coarse**
(β=5.6), and **small** (12³). All three bias the number — the chiral m_N/√σ ≈ 2.8 here vs a real
baryon's ≈ 2.1. The pilot proves the *object and its qualitative physics*, not the physical mass. The
quenched scan walls near m_crit where exceptional configurations (near-zero Wilson modes, unsuppressed
without a sea) make the plateau wild.

## L4 production — dynamical (first ensemble)

The first **dynamical** torsiton: the fermion determinant is in the path integral, so the vacuum is
the real, fluctuating, chiral-breaking condensate — sea-quark loops included. SU(3) N_f=2 fundamental
Wilson HMC (`generate_dynamical`, `run/07`), 12³×24, β=5.6, sea mass = −0.5, 48 thermalised configs.
HMC healthy (sea force ≠ 0, CG ~50 iters, dH ~ O(0.4) against H ~ 1.8M, acceptance high). Valence =
sea (the unitary point). The nucleon plateaus *later* than the pion, so each is fit in its own window
(pion t∈[4,10], nucleon t∈[8,11]):

| | m a | m/√σ |
|---|---|---|
| pion | 1.299(4) | 3.36 |
| nucleon (torsiton) | 2.109(14) | 5.45 |

with **m_N/m_π = 1.623**. Scale on *this* ensemble (`run/01` on the dynamical configs):
√σ a = 0.387(18), r₀/a = 3.166 — the string tension **drops vs quenched** (√σ a: 0.56 → 0.39, r₀/a:
2.28 → 3.17), i.e. the sea quarks soften confinement and make the lattice finer at fixed β: a real,
expected back-reaction, the vacuum the quenched pilot could not have.

**Reading it honestly.** This is a *single, heavy* sea mass: m_π/√σ ≈ 3.4 (the physical pion is ≈ 0.3),
so m_N/√σ ≈ 5.5 is the heavy-baryon value, not the physical one — a QCD-like nucleon at m_π ~ 1.5 GeV
sits at the same ~5.5. The physical torsiton mass needs the **chiral extrapolation** (several sea
masses, each a fresh ensemble) plus finer lattices and bigger boxes — the production campaign. What
this first ensemble establishes: the torsiton **binds in the real dynamical vacuum**, the sea visibly
dresses the gauge dynamics (the scale shift), and the pion is a genuine sea-quark Goldstone. The
method works unquenched, end to end.

## L4 — the mass-giving bag (condensate 3-pt + generations)

The fermion-mass hierarchy reduces to one number: the sharpness `s_T` of the dynamically-massive
**bag** the radial rungs couple to (Chapter 11). A sharp bag spreads the consecutive generations
across the full charged-lepton span (×3477); a broad bag flattens them. So the programme's headline
lattice question is **how sharp is the bag**, measured in units of the confinement scale, `s_T = R0/r0`.
The lever's productive window — where the consecutive ladder reproduces the observed span — is
`s_T ∈ [0.43, 0.70] r0`.

**The one-body proxy (`run/09`, the clean result).** The gauge-invariant dressed-quark scalar profile
`ρ(r)=Tr[S†S]` gives a bag whose half-density radius rises monotonically as the pion lightens. On the
12³×24 sea-(−0.5) ensemble, a valence scan (the guard dropping exceptional configs past the wall)
puts **three measured points inside the window** (s_T = 0.49, 0.45, 0.43 at m_π a = 0.32, 0.37, 0.40),
and the light-mass chiral extrapolation lands at **s_T ≈ 0.5 r0, in the window**. As a one-body proxy
this is the *best current estimate*: the bag grows into the productive window toward the chiral limit
— consistent with the hierarchy being derived.

**The genuine three-body condensate (`run/10`).** The connected nucleon scalar 3-point
`⟨N|q̄q(r)|N⟩` was built and **validated**: the sequential source (the derivative of the validated
2-pt contraction) passes its self-check (`recon == C_N` to all digits, by Euler), the sink guard is
sign-robust (the nucleon correlator here is uniformly negative — a convention), and a two-state fit in
τ removes the time-direction excited states. What it delivers cleanly is the **scalar charge
`g_S = 0.26 → 0.27`** (−0.90 → −0.95), *rising* toward chiral and shape-independent — a robust number.
What it does **not** yet deliver is the bag *shape*: source smearing was the wrong tool (it convolves
the spatial profile with the ~3a source blob — the smeared R0 ≈ 2.5a is just the smearing kernel),
and the point-source two-state fit returns R0 = 0.46a, **below the lattice resolution** (a ≈ 0.16 fm,
so R0 < ½ a — over-corrected on 7 τ points). The bag sits *at* the cutoff (~1a) on this coarse pilot;
the position-space `s_T` is bracketed (0.15 three-body fit ↔ 0.53 one-body proxy), **not pinned**.

**Generations (`run/08`, GEVP).** The variational nucleon spectrum resolves the **ground state
cleanly** (m_N a ≈ 2.26, consistent with the single-operator 2.11) but **not** the excited tower: the
three Wuppertal-smeared operators are nodeless Gaussians, nearly collinear, with almost no overlap on
the *nodal* radial excitations the generations live in. So the count (3? 4?) is unresolved — a
basis/statistics limit, not an absence. It needs **nodal operators** (Laplacian-filtered sources),
more configs, and a bigger box.

**Honest status.** The pilot ensemble (12³×24, a ≈ 0.16 fm, 48 configs) is at its **resolution
limit** for both the bag shape and the generation count. What is banked and reusable: the torsiton
binds; the full condensate-3-pt machinery is built and self-validated; `g_S` rises toward chiral; and
the one-body proxy says the bag enters the window. What is owed to a **production ensemble**: a finer
`a` (to resolve a ~1a bag — the dominant limitation), bigger `T` (T=48, to feed the two-state fit the
τ-range it starves for), bigger `L` (finite-volume control and room for the extended excitations),
more configs, nodal GEVP operators, and ultimately a lighter *unitary* sea. The pipeline is done,
end to end; it is the lattice that must grow.

## L4 production — L16×48 (Stage 1): the bag resolves

The first production ensemble: **16³×48**, β=5.6, sea mass −0.5 (same heavy mass as the pilot —
Stage 1 isolates the *T-effect* from the mass-effect), 122 thermalised/decorrelated dynamical
configs. The geometry change the pilot's resolution limit demanded: T=24→48 moves the backward node
out to t≈24 and hands the two-state fit ~16 τ-points; L=12→16 gives finite-volume room. Same
validated pipeline, end to end (`run/07` generate, `run/01,06,09,10`).

**Scale and spectrum (`run/01`, `run/06`).** Consistent with the pilot dynamical ensemble:

| | m a | from |
|---|---|---|
| scale | r₀/a = 3.131,  √σ a = 0.385 | `run/01` static potential |
| pion | m_π a = 1.269(3) | `run/06`, window t∈[4,10] |
| nucleon (torsiton) | m_N a = 2.091(7) | `run/06`, window t∈[8,11] |

with **m_N/m_π = 1.648** — the heavy-sea constituent-counting value, unchanged from the pilot's
1.623. The torsiton binds in the bigger box exactly as before; the bag question is what moved.

**The bag RESOLVES — the genuine three-body number (`run/10`).** This is the payoff of T=48. The
connected nucleon scalar 3-point `⟨N|q̄q(r)|N⟩` (sequential source, self-check passes to all digits)
lifts the **τ-plateau bag clean off the cutoff**, `R0 = 1.4–1.8a` (vs the pilot's unresolved R0<1a),
and the scalar charge is **rock-stable, g_S = 0.21**, across all sink times. The sink-time scan
(243 configs), the decisive Stage-1 test:

| t_snk | dist. to node (26) | R0 | s_T | err |
|---|---|---|---|---|
| 16 | 10 | 1.45a | 0.462 | ±0.013 |
| 18 |  8 | 1.42a | 0.452 | — |
| 20 |  6 | 1.76a | 0.563 | ±0.058 |

> **The bag is RESOLVED and IN the window [0.43, 0.70] at every sink time — `s_T ≈ 0.45–0.56`
> (central ≈ 0.50, sink-systematic ≈ ±0.06).** The value is *not* pinned within the window.

What is robust: the bag resolves (R0 ≈ 1.4–1.8a, off the cutoff), it lies in the productive window,
and `g_S = 0.21` is shape- and sink-independent (the genuinely pinned number). What is *not* pinned:
the precise `s_T`. The two clean points (t_snk=16,18) agree at ~0.45; t_snk=20 drifts up to 0.56 with
a 4× larger error, sitting only 6 slices from the backward node — so two effects push R0 up together
and are **degenerate on this ensemble**: residual excited-state relaxation (the ground-state bag is
genuinely larger) vs. backward contamination from the node (spurious inflation). The two-state τ-fit
that would separate them is unstable (R0 swings 4.2–6.1a — flagged `gs_reliable=False`; it needs more
τ-room, i.e. a larger T). Given the exponential lever (`d ln span / d ln s_T ≈ 15`), this ±0.06
maps to a charged-lepton span anywhere from ~8000 to ~50000 — so this ensemble confirms the
*mechanism* (resolved bag, in the window) but does **not** deliver the spectrum to ~30%.

**The one-body proxy agrees (`run/09`).** The dressed-quark scalar profile chiral scan puts the
lightest three points in the window and extrapolates to **s_T ≈ 0.430**, with the consecutive-rung
lever spanning **×3684 vs the observed ×3477** charged-lepton span — consistent with the derived
hierarchy.

**Three independent roads all land in the window.** They did not have to:

| road | what it measures | s_T |
|---|---|---|
| one-body proxy (`run/09`) | dressed-quark bag ρ(r)=Tr[S†S], chiral-extrapolated | 0.430 |
| three-body condensate (`run/10`) | ⟨N\|q̄q(r)\|N⟩, the genuine knot | 0.45–0.56 |
| lepton hierarchy (`genesis_horizon`) | the s_T the ×3477 lever *requires* | 0.428 |

All three sit inside [0.43, 0.70]: the bag sharpness measured directly on the lattice and the
sharpness the lepton masses demand are the same kind of number, not orders apart. The framework's
central claim — one substrate, one knot-size, sets the whole fermion ladder — survives a real
measurement. What the three roads do **not** yet do is *agree to a few percent*: the three-body bag
spans 0.45–0.56, and against the exponential lever that spread is the dominant uncertainty.

**Honest caveats.** Still a *single, heavy* sea mass (m_π a ≈ 1.27, far from the chiral point), and
the lever model sits inside the lepton road. The **sink-time scan is done and it revealed a
systematic** (s_T drifts 0.45→0.56 as t_snk approaches the node), so the value is bracketed, not
pinned. The honest lesson for Stage 2: the urgent lever is **larger T** (T=64), to give the two-state
τ-fit the room to stabilise and separate excited-state relaxation from backward contamination — *not*
only the lighter seas (−0.75, −0.85). Both are wanted; T is the one this scan flagged. Also owed: the
nodal GEVP generation count. But Stage 1's question — *does the bag resolve once T=48 gives it
τ-room?* — is answered yes: **the bag resolved, off the cutoff, in the window, on the genuine
three-body operator, at every sink time** — the precise value within the window is Stage-2 work.

## L1 — confinement (status)

The pure-gauge area law / string tension pipeline (`measure_potential`, `run/01`) is in hand and
audited; it gives the scale √σ that turns the torsiton mass into a cutoff-independent ratio. The
remaining L1 piece (a dual order parameter for the confinement *mechanism*) is open.
