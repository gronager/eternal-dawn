# Lab note — the information bounce (bring in C: entanglement) (2026-06-14, morning)

Status: **lab notes, brainstorming.** Forward-looking companion to
`spin-sector-and-the-bounce-source.md` (which retired "matter is torsion"). This one follows the
thread that *replaced the foundation*: the bounce as an information/entanglement event.

## The move
A (no-singularities + bounce → universes) + B (holographic/information ceiling) gave a *limit* but no
*force*. C = **quantum entanglement + unitarity** supplies the force. Two directions, both grounded:
1. **"Information cannot die" = unitarity.** Upgrades the no-singularities *axiom* to a *theorem*: a
   singularity destroys information, unitarity forbids it, ∴ a bounce.
2. **"The bounce is information fighting back" = the BH information paradox, resolved by the bounce.**
   Collapsing information, forbidden to be destroyed or packed past A/4, bounces into a child universe.

Grounding (not hand-waving): Bekenstein–Hawking S = A/4 *is* entanglement entropy; ER=EPR /
"entanglement is the glue of geometry" (Van Raamsdonk, Maldacena–Susskind); Jacobson's entanglement
*equilibrium* → Einstein's equations; baby universes / islands as an information-paradox resolution.

## The decisive computable result — the void seed is pure information
The void seed falls out of "matter's bit-count = holographic bound," **with no torsion**:
N fermions saturate A/4 = πR_s²/l_P² at N ~ (M_Pl/m_n)²/4π ~ 1.3e37, M_seed = N m_n ~ M_Pl²/m_n
~ a mountain mass. So the void seed comes out **three independent ways**:

| route | uses | M_seed |
|---|---|---|
| torsion bounce ρ_C | G, m_n, ℏ | M_Pl²/m_n |
| holographic/Schwarzschild | G, m_n | M_Pl²/m_n |
| **information saturation** | **S = A/4, no torsion** | **M_Pl²/m_n** |

One number, three derivations, the simplest needs no torsion ⇒ **the bounce scale is informational;
torsion was scaffolding.** `sims/.../holographic_seed.py` carries the first two; the third is the
N=(M_Pl/m_n)² one-liner above.

## The candidate unification
- **Bounce** = information resisting *destruction* (can't be over-packed past A/4).
- **Particle** = the dual: information refusing to *disperse* — a stable self-bound entanglement knot
  (a soliton / a bound level in a well). One principle (quantum information), both jobs — where
  torsion stranded into two sectors.

## The frontier calc (the whole ballgame): does the bound push back?
Need: the **sign of the back-reaction** as a collapsing ball's entanglement entropy crosses A/4.

Program:
- Track three curves vs R (fixed N): S_matter(R)~N (volume, ~fixed); S_bound(R)=A/4=πR²/l_P²
  (area, ∝R²); they cross at R_sat = the void seed.
- Build F(R) = E(R) − T_ent(R)·S(R), S capped at S_bound, T_ent ~ ℏc/(2π k_B R) (modular/entanglement
  temperature). Pressure P = −dF/dV. Back-reaction sign = sign of P (or V_eff'') near R_sat.
- **Tier-1 first pass (entanglement first law).** As R shrinks past R_sat the bound capacity falls,
  so the matter must shed dS_bound = (2πR/l_P²)dR; the energy cost dE = T_ent dS_bound = (ℏc/l_P²)dR;
  pressure P_ent = dE/dV = ℏc/(4π R² l_P²). At the void seed this is **REPULSIVE and ~ the matter
  energy density** (P_ent/ε ~ O(1–10)) ⇒ an information-pressure bounce is *plausible*, not excluded.
  (O(1) factors and the ρ_C coefficient ~1e7 are unpinned; sign is the robust part.)

Honest difficulty: saturating A/4 = forming a black hole = where classical GR makes a singularity.
A *bounce* there is intrinsically a quantum-gravity/information statement. The rigorous version is the
**generalized entropy** S_gen = A/4 + S_matter and the **quantum extremal surface** (Page-curve /
island machinery, 2019+): the bounce = the QES preventing the singularity, with information split
between the child universe (island) and Hawking radiation (the framework's CMB). And a real
sub-question: **bounce vs. radiate** — does the collapse hit the bound faster than entropy can Hawking-
radiate away? Fast ⇒ trapped ⇒ bounce; slow ⇒ ordinary evaporation. (Framework wants both: child
universe inside + Hawking CMB outside = the Page curve.)

Tiers: (1) toy thermodynamic V_eff / entropic pressure [sign, doable now — first pass above];
(2) Jacobson entanglement-equilibrium second variation; (3) generalized entropy / QES (rigorous).

## Results (Tier 2 + Tier 3 + the other limit) — added 2026-06-14

**Tier 2 (the rigorous sign).** Jacobson's first variation = Einstein (attractive). The bounce force is
the *entropy bound*, now a theorem-grade inequality: Casini (2008) proved the Bekenstein bound from
**positivity of relative entropy** (S_rel = ΔK − ΔS ≥ 0 — "information cannot die," as an inequality);
Bousso + Wall (GSL) give the holographic A/4 form. The bound is **one-sided**, so over-compression past
saturation is *forbidden*, not softly resisted → the back-reaction sign is **definite: repulsive wall.**
Chain: unitarity → relative entropy ≥ 0 → entropy bound → no compression past the void seed → bounce.

**Tier 3 (the QES = the child universe).** S_gen(r) = πr²/ℓ_P² + S_matter(r). The area capacity shrinks
(∝r²) as the collapse proceeds while S_matter ~ N stays; the quantum extremal surface sits where they
meet, R_QES = √(N/π) ℓ_P with N = (M_Pl/m_n)². Numbers: **R_QES = 0.12 fm ≈ void-seed R_s (0.42 fm) ≈
Compton (0.21 fm) — same scale.** The island caps the entropy at A_QES/4G = N (finite, holographic) =
the void-seed BH entropy — instead of running away. So **the island = the child universe**, the **QES =
the bounce membrane**, the **bounce = the Page-curve turnover** (info preserved), and **black hole =
child universe (island) + Hawking radiation (the CMB)** — the framework's central claim, *verbatim*, in
the language that resolved the information paradox. (Heuristic transcription of the island formula to a
collapsing ball in our universe; O(1)/coefficient factors unpinned; rigorous version = the gravitational
path integral. Scale and logic are sound.)

**The other limit — fermion birth (the dual).** `curved_gap` (M = √(M_flat² − R/4)) *is* the information
statement: the chiral condensate is an entanglement order parameter; curvature (Lichnerowicz R/4)
disrupts it, low curvature lets it crystallize. So the two limits are one knob (curvature), two ordered
phases: **bounce** = information forbidden to be *destroyed* (holographic ceiling → wall); **particles**
= information forbidden to *disperse* (crystallizes into the condensate, then the bound-level spectrum).
The bounce membrane R_crit = 4M_flat² is the phase boundary; it's persistent (the ongoing dawn) because
the parent holds the inflow on the high-curvature/fluid side and it crystallizes as it crosses.

**The arc.** "Is the matter torsion?" → torsion strands into two sectors and can't make the masses
without a postulate → the bounce, the void seed, the no-singularity axiom, the force, and the child
universe all re-derive from **quantum information / entanglement**, with torsion demoted to one possible
microscopic skin on the gravitational sector. Conserved throughout: the information never died.

## Bounce dynamics — the soft wall is the Planck force (added 2026-06-14)
No hard bounces: a hard (infinite) wall is a *discontinuity* (forbidden by the no-discontinuity axiom),
AND an inelastic bounce *destroys information* (forbidden by unitarity). So the bounce is **elastic**.
The Tier-1 resistance dE = (ℏc/ℓ_P²)|dR| gives a restoring force **F = ℏc/ℓ_P² = c⁴/G — the PLANCK
FORCE** (the GR maximum-force bound, Gibbons/Schiller). Finite ⇒ elastic, no discontinuity; maximal ⇒
the hardest elastic bounce nature permits. Max compression past saturation: ΔR = E/F = GM/c² = R_s/2
(~half the Schwarzschild radius — a real golf-ball squish, 0.21 fm for the void seed). The bounce is as
stiff as physics allows *without* a discontinuity: the **supremum-force bounce**. Both core axioms (no
discontinuities; unitarity/ongoing dawn) say the same thing — the bounce is soft and elastic.

## The hierarchy, and what's universal vs local (added 2026-06-14)
- **m_n/M_Pl = exp(−44).** NOT tuned (exp of an O(40) number). The **44 = ln(M_Pl/m_n) = e-folds of
  cohesion running = ½·ln(void-seed bits) = 1/g_UV** (up to the β-coefficient). **Located but underived:**
  "information cannot disperse" gives the *structure* (stable knot ⇒ soliton ⇒ dimensional transmutation,
  scale = M_Pl·e^(−1/g)), not the *number*. The 44 is an **input**. Sanity check: Diósi–Penrose
  gravitational self-localization sits at M_Pl, not m_n — so gravity doesn't set the matter scale either;
  cohesion does.
- **Why that 44 = selection, not derivation:** Smolin's CNS (already in the lineage) tunes it for
  **fecundity** = max black-hole/bounce production = **max information-reproduction**. An evolutionary
  *why*, honestly labeled.
- **Two kinds of "constant":** UNIVERSAL & DISCRETE (topology — charges, spins, gauge structure, particle
  *types*, and the generation *count* as a function of well depth) vs LOCAL & CONTINUOUS (scales — masses,
  set by the 44). **Mass = scale = flows; charge/spin = topology = fixed.** An electron is an electron
  everywhere (charge −1, spin ½) but weighs what its universe's 44 dictates.
- **Handshake / chemistry corollary:** matter carried into a sibling universe rebinds to *that* universe's
  local fields (the constants are local, not carried with you) → atoms resize, bonds re-energize → you
  disintegrate trying to shake hands. BUT note the nuance: if the dimensionless **ratios** (m_p/m_e, the
  12-spectrum) are **canonic** (universal well shape), chemistry is the *same structure, rescaled* (same
  periodic table, different absolute sizes/energies); only if the ratios are **selected** (well shape
  varies) is chemistry *qualitatively* different. Open which.

## Eq 6.2 in the age of information — scale vs spectrum (added 2026-06-14)
- **The scale (the 44) is measured, not derived** — one input. The weltformel is fixed up to *one* mass.
- **Are the ratios (m_p/m_e=1836, all 12) canonic or local?** This reduces to: **is the well *shape*
  (s_T, the coupling ratios) universal?** The framework's natural answer is *yes* — the well is the
  **self-consistent soliton** of eq 6.2 (the field digs its own well via the condensate); the universal
  field equation ⇒ a universal *dimensionless* shape ⇒ canonic ratios. Then: **measure 1 mass (scale) →
  predict all 12.** The coupling *ratios* (λ, g_v, the Fierz G_V/G_S=2) are already scheme-independent
  (`coupling_derivation`), supporting canonic.
- **Strong-sector ratios ARE canonic** (QCD-like dynamics, lattice-predicted — same in any such universe).
- **The unsolved crux:** the cross-sector m_p/m_e and the full generation hierarchy (~3477× lepton span)
  are NOT cleanly reproduced by the well (`dirac_woods_saxon`: ~O(1–10), not 3477). So "canonic spectrum
  from one universal well" is the **goal**, and the generation mechanism that would make it true is the
  open frontier. A Woods–Saxon FIT (scale + one ratio → 12) is descriptive; the *physics* that fixes the
  shape canonically is owed. The 12 = **4 towers (charge/colour grip, c_T) × 3 rungs (radial levels, s_T)**;
  s_T sets the level spacing, c_T the tower split — both should be canonic if the soliton shape is.

## The condensate is primary; torsion is a UV skin; the vector channel is the grip (added 2026-06-14)
The eq-6.2 audit + the "is it still scaled spin/torsion?" question resolve together:
- **The binding *is* the chiral condensate σ=⟨ψ̄ψ⟩** — one object, two names. "Spin/torsion" is a claim about
  its UV *origin* (dead as gravitational torsion: 37 orders subcritical → it's strong IR dynamics by
  dimensional transmutation, whatever the deep origin). "Information-dispersion resistance" is what σ *is*
  in principle: an entanglement order parameter / off-diagonal long-range order (the BCS analogy the
  chapter already makes — Nambu's "mass is a gap" = "mass is sustained coherence"). Not competing operators;
  a UV story and a principle story for the same σ. The scalar well S=−gσ is unchanged either way.
- **What the information view changes: what's *primary*.** Torsion forces the axial operator (ψ̄γ₅γ^μψ)²
  as fundamental → the scalar well is a derived Fierz shadow (dragging the mandatory vector partner,
  G_V/G_S=2). Information inverts it: the scalar condensate is primary, and "won't disperse" is a
  **gradient** statement (∇σ)² — the bag's **surface tension**, already in the model (Friedberg–Lee σ-field,
  `cs._sigma_newton`). ⇒ the contact term eq:nld is the integrated-out UV artifact; the σ-field-with-
  stiffness (eq:gap onward) is the load-bearing IR description, and the "cohesion" *is* the surface tension.
- **The eq-6.2 missing term (now fixed): the vector repulsion.** eq:nld is axial-axial; the bridge to the
  scalar well is a Fierz step that also yields a vector (ω-like) channel, G_V/G_S=2 (`coupling_derivation`,
  scheme-independent, from m_σ=2M, m_ω=√6 M). The book's eq:radial kept only the scalar S(r); the vector
  V(r) (entering E−V₀) was dropped. It's the term that *unbound* the soliton in the 2nd-order reduction;
  binding survives only in the first-order (G,F) solver (`dirac_soliton`), where the vector **loosens the
  bag monotonically**.
- **The vector channel is the m_e/m_p freezing knob.** g_vω is the conserved-charge (baryon/vector)
  coupling: a baryon feels scalar+vector, a colourless lepton doesn't feel the strong vector → the
  scalar/vector balance *is* the tower-splitting grip c_T. Strong-sector internal ratios are frozen by
  dimensional transmutation (one Λ → all hadron ratios universal — not a hope, how QCD works). Cross-sector
  m_e/m_p is frozen **iff** the grip is universal — and G_V/G_S=2 is scheme-independent (a pure number),
  which is the evidence that it is. So **m_e/m_p frozen ⇔ G_V/G_S a topological number** = the lab-note
  "topology fixes ratios, scale floats," now wearing an operator (vector=charge/topology side, scalar=
  mass/scale side). Freezing is argued; the *value* 1836 is not yet reproduced (mean field gives O(1–10)).
  If the lattice confirms G_V/G_S≈2, that is direct evidence m_e/m_p is the *same* in a sibling universe.
- **Numerical check (`dirac_soliton`, adiabatic g_v ramp).** Converged window g_v≤0.5: every level shifts
  up, spacing ratio d₂₃/d₁₂ compresses 6.66→6.38 — the spectrum is *not* the scalar well's alone. Past
  g_v≈0.75 the **third bound level is ejected** and convergence breaks — the generation *count* rides on
  the scalar/vector balance, not just the core sharpness. (The fragility is the falsifiable content.)

## Measure the well, don't parametrize it — count vs. span, and the wall is the lever (added 2026-06-14)
Built the shape-agnostic well→spectrum pipeline (`well_spectrum.py`): consumes any (M(r),V(r)) — mean-field,
parametric, or a lattice-measured profile (`well_from_bag_profile` maps the dressed-quark density ρ(r)
to M(r)=m_vac(1−ρ)) — and returns the spectrum as *outputs* (level count, node-labelled E_n, two
configmass pieces). Engine: inward/outward RK4 matching with a **pole-free Wronskian** eigen-condition.
- **Doubler finding (matters for the chapter).** The old finite-difference matrix solver was
  fermion-doubler-contaminated: grid-scale spurious modes were counted as radial excitations and
  **inflated the generation count**. The doubler-free count is smaller; any "the bag binds 3" resting on
  the matrix count was unreliable. The un-fitted mean-field well binds **1** level at standard couplings.
- **The slip, and its correction.** Went "don't parametrize with 2-param Woods–Saxon" → "measure the full
  M(r),V(r)" → (regressed) "measure two numbers s_T and m_vac·r0." That last step is the *same* 2-parameter
  description, no better than WS. **Test (fixed R0 & depth, 5 different shapes):** the **count** is robust
  (11–12 everywhere — a coarse WKB-area, ~2-number feature), but the **span swings 23×** (3.4 → 79). So:
  count ≈ 2 numbers (fine); **masses/span = full shape (2 numbers far too light).**
- **The wall/tail is the lever, not the radius.** At fixed R0 the span is set by the surface thickness:
  sharp wall (a=0.05R0)→span 79, soft (a=0.45R0)→3.4. The chapter's s_T=R0/r0 *ties* a=0.15R0 and so
  conflates size with wall-sharpness; the honest lever is **how fast M(r) falls at the surface**. The
  lattice must **resolve the wall/tail**, not just the half-density radius.
- **The two measured numbers that ARE legitimate, and what they fix:** m_vac·r0 (depth) fixes the **count**
  (at fixed shape: 1→2→3→7→12 as m_vac·r0 = 0.9→2→4→8→16) and the overall scale; s_T (size) is coarse.
  The **hierarchy** is carried by the wall, which no two numbers capture.
- **Derived-not-fitted = the line between light and real.** What makes "measure the full M(r) and solve"
  NOT just WS-with-more-knobs is *where M(r) comes from*: the **dressed-quark self-energy output by the
  gauge dynamics at the one scale Λ**, measured from the *propagator* **independently of the lepton masses**,
  then used to *predict* them — the same status as lattice QCD predicting hadrons from Λ_QCD. Light = fit
  the shape to the answer; not-light = compute the propagator, then solve. **Open hard part (not papered
  over):** the closed-form M(r) from Λ alone exists only at mean-field, where it's too shallow (binds 1);
  the lattice propagator is the non-perturbative *evaluation*, not a closed-form derivation.

## L16×64@−0.75 clean depth: the bag binds ~one quark, not a tower (added 2026-06-14)
First clean measurement on the live ensemble (89 configs, NO_MPS=1 direct-GPU; the MPS client-connect
bug needed the empty-pipe bypass): **m_π a = 0.6336(20)**, **m_N a = 1.1159(79)**, m_N/m_π = 1.76. So:
- **Depth m_vac·r0 = (m_N/3)(r0/a) ≈ 1.16** (r0/a=3.131 from L16×48, to be pinned by run/03). With the
  L16×48 bag shape this binds **0–1 levels** — not three.
- **The robust statement (WKB):** n_levels ~ m_vac·R0/π, with m_vac·R0 = (m_vac·r0)·s_T ≈ 1.16·0.6 ≈ **0.7**,
  while 3 levels need m_vac·R0 ≳ 2.5π ≈ **8**. The measured QCD-like bag is **~11× short** of a three-rung
  tower. Physically sensible: a dressed/constituent quark is a *ground-state* object; its radial
  excitations are unstable resonances (ρ′, ρ″), not three stable rungs.
- **So the "3 generations = 3 radial levels in one bag" picture is in ~order-of-magnitude tension with the
  lattice.** Caveats (keep it a tension, not a verdict): heavy pion (m_π≈780 MeV, not chiral); the −0.75
  bag *shape* not yet measured (used −0.5); r0/a and the M=m_vac(1−ρ) dictionary assumed. But an 11× gap
  in m_vac·R0 is not closed by modest chiral trends. **This is the sharpest quantitative test yet of the
  generation-tower picture, and it is failing.** → revisit Picture B (distinct solitons per generation)
  and the decay-diagram structure (separate conserved lepton flavors ⇒ generations are not radial
  excitations of one well — see next).

## Generations are topology, not excitation — the decay diagrams settle it (added 2026-06-14)
Independent of the lattice, the μ/τ decay diagrams kill the radial-tower picture and point to distinct
topologies. Concordant with the lattice (m_vac·R0~0.7, no tower), now from particle data.
- **The killer: μ→eγ does not happen (BR < 4.2e-13).** If μ and e were radial levels (2s,1s) of one
  well, the muon would de-excite radiatively to the electron — same J=½, same parity, an M1 transition,
  *fast*. It refuses; it only dies *weakly* (μ→ν_μ W, W→eν̄_e). Radial excitations carry no protected
  label; the absence means the generation label is a **conserved quantum number (lepton flavor)**, so the
  generations are **NOT radial excitations of one well.**
- **The vertices:** μ→ν_μ is within-generation (the lepton converts to its *own* neutrino); the W carries
  the flavor change into a *separate* vertex. e and μ never touch directly. τ decays *democratically* to
  e/μ/hadrons (≈18/17/65%) — no rung-by-rung ladder (another nail).
- **Reframing (the framework's own "knot" language): generations = distinct knot TOPOLOGIES.**
  Topologically protected ⇒ can't smoothly deform into each other ⇒ no μ→eγ; flavor conservation is
  *automatic*, not imposed. The **W is the only knot-reconnection operator** (cuts/re-ties the topology,
  carrying flavor away); strong/condensate dynamics binds the knot, only the weak interaction re-ties it.
- **Two topological labels:** charge (shared by all 3 generations, sets the universal W-coupling) ⊥
  generation (distinct, sets the mass via knot complexity). Extends "charge=topology=fixed, mass=scale=
  flows" with a *second* fixed index. Universality (W-coupling identical across gens, <1%) = the
  generation label is orthogonal to charge.
- **Neutrino prediction (the handle):** charged knots are pinned (no mixing, no μ→eγ); neutral knots
  (neutrinos) are light and near-degenerate ⇒ **tunnel between topologies ⇒ mix** (oscillations). Predicts
  the KNOWN anti-correlation the SM merely inputs: **mixing ∝ 1/hierarchy** — charged leptons steep
  (3477) → ~no mixing; quarks medium → small mixing (CKM); neutrinos flat (near-degenerate) → large
  mixing (PMNS). Same structure that forbids μ→eγ permits ν-oscillation.
- **Majorana/Dirac:** a neutral knot is its own antiknot ⇒ **neutrinos Majorana**; charged knots ≠ their
  antiknots ⇒ **charged leptons Dirac**. Falsifiable via neutrinoless double beta decay.
- **Why this helps:** the hierarchy is no longer "level spacing in one well" (lattice-failed, capped at
  span ~3–5) but the **field energies of distinct soliton topologies** — a more complex knot can be far
  heavier, *unbounded* by any single-well level spacing. Doesn't hand us 3477, but removes the ceiling.
- **Lineage (resembles, not adopts):** Bilson-Thompson braided preons, Harari–Shupe rishons — generations
  as topology, with the decay diagrams as the same evidence.

## Start with the neutrino; W/Z/Higgs are waves not knots (added 2026-06-14)
- **Neutrino = the bare generation topology.** Matter = a topological soliton of the condensate is already
  real physics: the **Skyrmion** (a π₃ winding of the chiral field, M=SU(2)≅S³, π₃(S³)=ℤ) *is* baryon
  number (Witten). Layer the charges as windings: baryon# = Skyrme π₃; electric charge = a gauged U(1)
  winding; **generation = a THIRD invariant** of the richer (chiral×axial/torsion) manifold. A charged
  lepton carries all three tangled; a **neutrino carries the generation invariant alone** (no charge, no
  colour, ~no mass) → it is the clean probe. *Start with the neutrino.*
- **Mixing ∝ 1/hierarchy is a PREDICTION (SM inputs it).** Mixing = near-degenerate knots tunnelling
  (2-level QM: angle ~ off-diagonal/splitting). Charged leptons steep (3477)→~0 mixing; quarks medium→
  small (CKM); neutrinos near-degenerate→large (PMNS). Same near-degeneracy that lets ν oscillate is what
  is absent for the steeply-split charged leptons (→ no μ→eγ). **Neutrinos Majorana** (neutral=own antiknot),
  **charged Dirac** — testable via 0νββ.
- **W/Z/Higgs are NOT solitons — the *other* excitation class.** A condensate supports knots (topological,
  conserved, stable = the fermions = stable information) AND modes (non-topological ripples = the bosons =
  radiated messengers). **Higgs = σ**, the scalar breathing mode of the binding condensate (sets every
  knot's mass; v→0 kills them; the m_σ~2M scale) — Peter Higgs keeps his particle, as the condensate's
  amplitude mode. **W/Z = vector modes + eaten Goldstones**; the **W is the knot-reconnection quantum** of
  the previous section (a wave that re-ties a knot, carrying off ΔQ and the flavour). Caveat: this is
  composite-EWSB → owes the **S-parameter** (S<0.1, the existing `electroweak_S.py` / L-targets), and the
  standing "NOT technicolor" caution (resembles the EWSB pattern, not the model-building).

## Homotopy first pass — what gives a 3-valued generation invariant on the neutral knot (added 2026-06-14)
The spine is solid (Skyrme: matter = π₃ texture, baryon# = winding). The generation needs a SECOND
invariant carried by the neutral knot. Three routes:
1. **Hopf / linking — π₃(S²)=ℤ.** If the neutral order parameter is an S² field n̂(x) (the chiral/axial
   direction), 3D textures are **Hopfions**: π₃(S²)=ℤ = self-**linking number** = "how knotted." Generation
   = linking number; the neutrino is the pure S² texture. Literally knots. Gives a *tower* (ℤ); "exactly 3"
   = energetics (only the lowest links stay bound — a finite knot tower). Mass ↑ with linking.
2. **Zero-mode index — Jackiw–Rebbi / Atiyah–Singer.** A knot binds an integer number of protected fermion
   zero modes = a topological index; generation = which mode. Could give **exactly 3** *if the index is 3*.
   The modes share the knot's charge (same charge across gens) but carry distinct conserved indices (flavour
   conserved → no μ→eγ). **Mechanism for mixing-vs-hierarchy:** the splitting of the zero-mode multiplet is
   sourced by the CHARGE coupling → charged knots split hugely (pinned), neutral knots stay near-degenerate
   (mix). This is the most attractive route (exactly-3 + the mixing correlation), *contingent* on the index
   = 3 — tying that to N_c=3 is a speculative NEW claim (N_c ⊥ N_gen in the SM; not established).
3. **Discrete flavour symmetry — π₀(ℤ₃ or A₄)=3.** Generation = which of 3 discrete vacua; connects to the
   established A₄/ℤ₃ neutrino-mixing literature (tri-bimaximal etc.); mixing = small symmetry breaking.
   Gives exactly 3 by construction, but a *postulated* discrete symmetry (less derived).
- **Honest verdict:** π₃ is generically ℤ (or ℤ^k), so a pure-topology "exactly 3" does NOT come from a
  finite homotopy group of the obvious manifolds — it comes from the **index** (route 2) or a **discrete
  symmetry** (route 3). Favoured synthesis: the knot is a Skyrmion (baryon# = π₃), its protected zero-mode
  multiplet is the generation (route 2), with the charge sourcing the splitting (→ neutral degenerate/mix,
  charged pinned). **The decisive calc:** count the fermion zero modes of the Eq-6.2 Dirac operator in a
  unit-winding (Skyrmion) background — is the protected multiplet 3-dimensional? (Honest: the B=1 Skyrmion's
  collective modes give the nucleon/Δ spin-isospin tower, NOT obviously 3 generations — so "index=3" is a
  hope to be checked, not a result.)

## Route 2 walked: one winding = ONE protected level, not three → pivot to Route 1 (added 2026-06-14)
Ran the existing CQSM instrument (`chiral_quark_soliton`: the Eq-6.2 Dirac in the hedgehog winding
background — the Route-2 zero-mode calc was already built; also `hedgehog`, `kahana_ripka`).
- **B=1 hedgehog → ONE topologically protected level** (K=0⁺, E/M=0.227 at physical depth) = the baryon
  (spectral flow: one zero-crossing per unit winding = baryon number). **Not three.**
- **The winding DOES help the tower bind (on depth):** excited K=0 levels appear; a 3-level window opens at
  M·r0~3 vs ~8 for the plain scalar well (~6× more efficient — the chiral rotation binds harder). Physical
  CQSM depth M·r0~1–2 → factor ~1.5–2 short (vs ~11× for the plain well). Supports the "winding/tower
  helps" instinct on the COUNT-depth front.
- **But the excited levels are NOT protected** — only the valence crosses zero (spectral flow); the rest
  are ordinary in-gap bound states = **Roper-like radial excitations** (N(1440)). In the real world the
  Roper DECAYS (N*→Nπ) — the strong-sector analog of the forbidden μ→eγ. So radial-tower states (winding or
  not) would decay; they can't be the stable generations.
- **Verdict:** Route 2 (zero-mode multiplet of ONE winding) gives 1 protected state, not 3. Protection
  (no μ→eγ) and tower (3 levels) are in TENSION: one winding protects one level; a stack of excited levels
  is unprotected (Roper-like).
- **Pivot to Route 1 (Hopf linking), which reconciles both instincts:** each linking number n (π₃(S²)=ℤ) is
  a DISTINCT topological sector → separately protected (no inter-generation transition, μ→eγ forbidden) AND
  an ordered ladder (mass↑ with linking = the tower). It's a tower *and* protected, because each rung is a
  different knot, not an excitation of one knot — exactly the distinction the Roper result forces. Neutrino
  carries link# alone → low links near-degenerate → mix. "Why 3" = energetics of bound links. **Next calc:
  Dirac zero modes in Hopfion backgrounds of link 0,1,2 — does the mass ladder up, and is each protected?**

## L16×64@−0.75 bag: scales up, doesn't sharpen — well picture not rescued (added 2026-06-14)
3-body condensate, 90 configs, t_snk=24, valence −0.5 on the −0.75 sea (g_S=0.234, self-check PASS).
- Measured bag: R0=2.71a, **s_T=0.855** (r0/a=3.165), wall **a_wall/R0=0.42**, tail λ=1.34a. Spectrum at
  the physical depth (m_vac·r0=1.16): **1 level**, span ~2.4 when forced. Same verdict as L16×48.
- **Wall trend L16×48→L16×64: a_wall/R0 0.46→0.42** (tail λ 0.96→1.34). The bag **scaled UP** (s_T 0.62→
  0.86) toward the lighter sea with the wall *fraction* ~constant — growth, not sharpening. The ~9% drop is
  marginal; both soft. **No rescue of the span.**
- Caveats (not final): valence −0.5 (not the −0.75 chiral point); the t_snk=24 plateau is flagged
  excited-contaminated (its unstable two-state fit reached for R0→0, so the true ground-state bag *could* be
  sharper — wants a larger t_snk / stabilised two-state fit). But as it stands the well/tower picture is not
  saved by the chiral trend → confirms the pivot to Hopfions (Route 1).
