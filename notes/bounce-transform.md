# The bounce transform: what actually happens at Cartasis

Working note. The question: when matter crosses the Cartasis membrane, what does
the geometry *do* to it — flip the arrow of time? flip charge/chirality (matter
↔ antimatter)? filter one handedness? Until this is pinned down, every
"inside vs outside" statement (CMB = Hawking, dark-matter filter fraction,
baryon asymmetry) rests on a verbal dictionary. The aim here is to replace the
dictionary with statements that are either theorems or clearly-labelled
postulates.

## 1. Three separable questions

The "bounce transform" conflates three logically independent things. Keep them
apart:

1. **Time orientation / arrow.** Does proper time reverse? Does the
   thermodynamic arrow reverse?
2. **Discrete charge symmetry (C).** Does matter emerge as antimatter?
3. **Chirality (handedness).** Is one handedness preferentially transmitted
   (filter), and can a net handedness be *created* vs merely *amplified*?

Our simulated bounce (Tier 1) already answers some of these; the rest reduce to
the C/P/T properties of the Einstein–Cartan interaction.

## 2. Arrow of time: a Janus point, not a reversal

In the Tier-1/OS bounce the matter worldlines carry a **single, monotonic proper
time** τ. The scale factor `a(τ)` has a minimum at the bounce; nothing runs
backwards in τ. So there is **no reversal of the geometric time orientation** —
the bounce is a smooth *continuation*, collapse → expansion.

The **thermodynamic** arrow is the interesting one. Track (gravitational)
entropy along τ:

- Parent side (τ<0): infalling matter came from a structured, high-entropy
  universe; compressing toward ρ_C dismantles structure and smooths the geometry
  (Weyl → 0). Entropy *decreases* toward the bounce.
- Baby side (τ>0): expansion lets structure (and eventually black holes) form
  again. Entropy *increases* away from the bounce.

So entropy has a **minimum at the membrane**, and the locally-defined arrow
(direction of increasing entropy) points *away from the bounce on both sides*.
This is a **Janus point** (Barbour–Koslowski–Mercati; cf. Aguirre–Gratton,
Carroll–Chen): one spacetime, monotonic proper time, two thermodynamic arrows
diverging from a central entropy minimum.

Consequence: the baby sees the bounce as its low-entropy past (a normal Big-Bang
arrow) — this is exactly the Past Hypothesis derived in Ch. 7, and it needs **no
time flip**. The "r becomes timelike / time inverts" language in Ch. 4 is a
coordinate statement inside the horizon, not a physical reversal of the matter's
proper time.

### The fork that matters

There are two geometrically distinct ways to read an entropy-minimum surface:

- **(A) Continuation (Janus).** The two sides are *physically distinct, real*
  regions (real parent, real baby) joined by continuous evolution. Time
  orientation is preserved; **C is preserved** (matter stays matter). This is
  what our dynamical bounce is.
- **(B) CPT mirror.** The surface is identified as a *CPT-reflection*: the "other
  side" is the CPT image of this side — reversed time orientation, parity-flipped,
  charge-conjugated → an **antimatter** partner universe. This is the
  Boyle–Finn–Turok "CPT-symmetric universe" construction (the Big Bang as a CPT
  mirror).

(A) and (B) are *different physics with different predictions*. The framework's
"Mechanism A: Flip (PT ⇒ C ⇒ antimatter)" is really scenario (B) — and it does
**not** follow from the Einstein–Cartan bounce dynamics. It would have to be
*added* as a discrete-symmetry postulate on the bounce surface. Our simulated
bounce is scenario (A).

## 3. C, P, T of the Einstein–Cartan interaction

This is the decisive part, and it is a symmetry theorem, not a model.

Integrating out torsion in Einstein–Cartan–Sciama–Kibble gives the Hehl–Datta
four-fermion contact interaction in the Dirac Lagrangian,

    L_HD = −(3/16) κ (ψ̄ γ⁵ γ_a ψ)(ψ̄ γ⁵ γ^a ψ),     κ = 8πG/c⁴,

i.e. minus the square of the **axial current** j⁵^a = ψ̄ γ⁵ γ^a ψ. Its discrete
symmetries:

- **C:** the axial current is C-**even** (unlike the vector current j^a, which is
  C-odd). So (j⁵)² is C-even.
- **P:** under parity j⁵^0 → −j⁵^0 and j⁵^i → +j⁵^i, so the contraction
  (j⁵_a j⁵^a) = (j⁵^0)² − (j⁵^i)² is P-**even**.
- **T:** likewise T-**even**.

Therefore **L_HD is separately C-, P-, and T-even.** The torsion interaction that
*drives the bounce* (the negative (axial)² term is exactly the −ρ²/ρ_C piece that
caps the density) carries **no intrinsic handedness and no charge preference.**

Three immediate corollaries:

1. **No flip.** A C-even interaction cannot turn matter into antimatter. Scenario
   (B) is not produced by EC dynamics; it is an extra postulate.
2. **No filter from nothing.** A C/P-even interaction cannot *generate* a chiral
   asymmetry. At most it can **amplify a pre-existing one**: in a background with
   net axial density ⟨j⁵⟩ ≠ 0, the mean-field interaction energy depends on
   alignment, so an existing bias is self-consistently enhanced — but ⟨j⁵⟩ = 0 in
   → 0 out.
3. **No arrow flip** (consistent with §2): T-even dynamics, continuous proper
   time.

So in **minimal ECSK the bounce is CPT-even**: it preserves the arrow (Janus,
not reversed), preserves charge, and cannot create chirality.

### Where a *derived* handedness could still come from

- **Parity-violating gravity:** add the Holst / Nieh–Yan term (finite
  Barbero–Immirzi parameter) or a non-minimal axial–torsion coupling. These are
  P-odd and *can* make torsion couple chirally. This is the honest place to look
  for a derived filter — and it is a specific, testable extension, not a fudge.
- **Standard-Model CP violation + sphalerons** at the trans-EW temperatures of
  the bounce: ordinary baryogenesis machinery operating in the bounce background.
- **Particle production at the bounce** (Poplawski): the bounce creates particles;
  the asymmetry rides on whatever CP violation is present, not on the EC term.

## 4. What this does to Mechanisms A / B / C (Ch. 4)

Current Ch. 4 "synthesis" = A (flip, dominant) + B (filter, modulation) + C
(statistical seed). The analysis above forces a revision:

- **Mechanism A (flip)** is **not** a consequence of the bounce. Demote it: it
  exists *only* if one adopts the CPT-mirror postulate (B-scenario,
  Boyle–Turok). Then it is a strong, separate hypothesis with its own
  predictions (e.g. an antimatter partner universe, specific neutrino content),
  to be evaluated on its own merits — not a derived bounce effect.
- **Mechanism B (filter)** is, in minimal ECSK, an **amplifier of a pre-existing
  bias, not a source**, and is symmetric (no preferred handedness) unless one
  adds parity-violating gravity. The "dark matter = filter-rejected antimatter,
  f ≈ 1/6" story needs either (i) a pre-existing bias to amplify, or (ii) the
  P-violating extension — state which.
- **Mechanism C (statistical seed)** is the **only** asymmetry *source* that
  survives in minimal ECSK, and it sits well with CPT symmetry of the supraverse
  (equal matter/antimatter lineages overall). This should be the load-bearing
  mechanism.

Net: the matter/antimatter asymmetry is **inherited (C)**, possibly **amplified
(B, mean-field)**, and the alternating-chirality "flip" picture requires the
explicit CPT-mirror postulate (A=B-scenario) rather than coming for free.

## 5. Why this must precede the CMB tests

The CMB = parent-Hawking identification and the "flip preserves thermality"
argument in Ch. 4 both lean on the transform. If the bounce is CPT-even
continuation (A), the CMB we see is the redshifted radiation of the *same-charge*
parent and there is no chirality distortion to look for; if it is a CPT mirror
(B), the partner radiation is C-conjugated and the predicted spectral/polarisation
signatures differ. We cannot interpret a CMB calculation as "great" or "bad"
without committing to A vs B first. **Settle the transform, then refine.**

## 6. Concrete next steps

1. **Janus entropy diagnostic** (cheap): add a coarse gravitational-entropy proxy
   to the Tier-1 bounce and show the minimum-at-membrane / diverging-arrows
   structure explicitly. (Homogeneous Weyl = 0, so use a perturbed or
   shear/structure proxy.)
2. **Hehl–Datta mean-field filter coefficient** (medium): solve the self-consistent
   ⟨j⁵⟩ amplification at ρ_C; quantify "amplifier not source," and get the
   amplification gain vs density.
3. **Parity-violating extension** (medium–hard): add the Holst/Nieh–Yan term and
   compute whether torsion then transmits one handedness preferentially — a
   *derived* filter, with a number to compare to f ≈ 1/6.
4. **Evaluate the CPT-mirror postulate** (theory): if the framework wants the
   flip, adopt Boyle–Turok explicitly and inherit its predictions; compare.
