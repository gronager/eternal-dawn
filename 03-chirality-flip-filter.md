# 03 — Chirality: Flip, Filter, or Both

## The question

When matter crosses Cartasis from parent to baby, what happens to its chirality (matter vs. antimatter content)? Three candidate mechanisms have been considered. The best current understanding is that two of them operate simultaneously.

## Mechanism A: Flip

**Claim:** The bounce geometry implements a CPT-like operation. Matter on the parent’s side emerges as antimatter on the baby’s side (or vice versa).

**Physical basis:** The bounce involves a smooth transition where the radial direction (timelike inside the parent’s horizon) becomes the time direction in the baby. This is a genuine geometric inversion of the time direction at the crossing. Combined with parity inversion that’s natural in spherical bounce geometry, this is a PT operation. By CPT symmetry (a theorem in any local Lorentz-invariant QFT), PT applied to matter is equivalent to C applied to matter — i.e., conversion to antimatter.

**Consequence:** Each generation has opposite chirality bias from its parent. Lineages alternate matter / antimatter / matter / antimatter across generations.

**Strengths:**

- Built on CPT symmetry, which is a theorem rather than an assumption
- Preserves total particle+antiparticle content, so every generation is viable
- Naturally produces matter-antimatter asymmetry without compounding to sterility
- Consistent with the observed near-thermal Planck spectrum of the CMB (assuming CMB = parent Hawking radiation, see `05-dark-sector.md`)

**Weaknesses:**

- The “genuine geometric inversion” claim is subtle; some readings of the bounce geometry don’t produce a physical chirality flip but only a coordinate relabeling
- Requires the bounce geometry to have specific topological properties that may not hold for all bounce solutions
- Hard to verify without detailed Einstein-Cartan bounce calculations

## Mechanism B: Filter

**Claim:** The torsion-spin coupling at Cartasis density is chirality-selective. One handedness experiences stronger effective repulsion than the other, given a background chirality bias. The minority chirality is preferentially scattered back; the majority chirality passes through.

**Physical basis:** The Einstein-Cartan coupling between torsion and fermion currents goes through the axial current $J^5 = \bar\psi \gamma^5 \gamma^\mu \psi$, which has definite handedness. At high spin density with a chirality bias, the effective interaction is asymmetric between handedness.

**Consequence:** Each generation has *enhanced* chirality bias relative to its parent. The filter compounds: small biases get amplified across generations.

**Strengths:**

- Cleanly grounded in the explicit form of the Einstein-Cartan Lagrangian
- Provides a mechanism for the strong asymmetry observed in our universe (deep-generation enhancement from many filterings)
- Can be calculated quantitatively from the Lagrangian

**Weaknesses:**

- If filter is the *only* mechanism, deep-generation universes approach pure-chirality limits where Standard Model physics breaks down (no antimatter for QFT loops, broken electroweak structure)
- Going up the lineage tree, parents have *less* asymmetry, eventually reaching near-symmetric universes that can’t form structure or black holes — the “hell upward” problem
- The lineage tree can’t get started from such an upper limit, breaking the framework’s self-consistency

## Mechanism C: Statistical inheritance from OG fluctuations

**Claim:** OG universes (those nucleated directly from vacuum fluctuations rather than descended from parents) inherit their chirality bias from the *specific* fluctuation that seeded them. The vacuum is symmetric on average, but individual large fluctuations have nonzero bias drawn from a statistical distribution.

**Physical basis:** A vacuum fluctuation that reaches Cartasis density and seeds a bounce is a specific quantum state with definite quantum numbers, not a vacuum average. The chirality bias of this state is drawn from the fluctuation distribution. By CPT symmetry of the vacuum, equal numbers of matter-biased and antimatter-biased fluctuations occur.

**Consequence:** OG universes come in matched matter/antimatter pairs (statistically). Each OG starts a lineage with its specific bias.

**Strengths:**

- No baryogenesis mechanism needed; the bias is statistical, not generated
- Preserves CPT symmetry at the supraverse level (equal matter/antimatter lineages)
- Provides initial conditions for the lineage dynamics

**Weaknesses:**

- Says nothing about how the bias evolves across generations within a lineage; needs A or B for that

## The synthesis: Flip + Filter + Statistical inheritance

The current best understanding is that all three mechanisms operate:

1. **OG bounces inherit asymmetry from vacuum fluctuations (Mechanism C).** Each OG starts with a specific bias.
1. **Subsequent bounces flip chirality (Mechanism A).** Lineages alternate matter/antimatter across generations. This preserves viability and avoids the hell-upward problem.
1. **Torsion coupling modulates the bias magnitude (Mechanism B).** Filter effects refine the bias slightly with each generation, but don’t dominate the dynamics.

Result: lineages look like alternating matter/antimatter signs with slowly varying magnitudes. Each universe in a lineage is viable. Going up the lineage tree, you find alternating-sign universes with comparable bias magnitudes, not progressively-symmetric hellish universes.

## The CMB constraint

The CMB is observed to be very nearly a perfect Planck spectrum (deviations < $10^{-5}$). If the CMB is interpreted as parent Hawking radiation crossing into our universe (see `05-dark-sector.md`), then a strong filter would produce specific spectral distortions (rejection of one chirality component during crossing). The observed near-perfect Planck spectrum constrains the filter to be weak.

A flip-dominated mechanism, by contrast, preserves the thermal spectrum (CPT conjugation doesn’t change photon thermality). So observations favor flip over strong filter.

This is consistent with the synthesis: flip is the dominant generational mechanism, filter is a small modulation.

## Implications for matter-antimatter asymmetry

In the synthesis:

- Our observed asymmetry ($\eta \sim 6 \times 10^{-10}$ baryon-to-photon ratio) comes primarily from the statistical bias of our OG ancestor’s seeding fluctuation, propagated through many generations of flip + filter modulation.
- The “missing antimatter” of standard cosmology is not actually missing — it exists in other lineages (matter-biased OG lineages and antimatter-biased OG lineages are equally numerous in the supraverse).
- Within our lineage, alternating generations have opposite bias signs. Our parent is antimatter-dominated. Our grandparent is matter-dominated (like us). Our children will be antimatter-dominated. Etc.

## Implications for dark matter

If filter operates even weakly, some matter content from the parent that doesn’t pass the filter remains at the Cartasis membrane. From the baby’s perspective, this is gravitational influence from “matter just outside our cosmological horizon, in the wrong chirality” — which manifests as dark matter (see `05-dark-sector.md`).

The dark-matter-to-baryon ratio (~5:1 by mass in our universe) might be quantitatively predicted by the filter efficiency. If our parent had bias $b_p$ and the filter passes fraction $f(b_p)$, then the dark matter fraction is $\sim (1-f)/f$. A 5:1 ratio implies $f \sim 1/6$, meaning the filter rejects about 5/6 of what crosses. This is a specific predictable number from Einstein-Cartan Lagrangian parameters.

## Open questions

- Precise calculation of flip vs. filter contributions from a specific bounce geometry
- Whether baryon number is exactly conserved at Cartasis, or violated through electroweak anomalies
- How quantum field theoretic effects (sphaleron processes at trans-EW temperatures) interact with the bounce dynamics
- Whether the chirality dynamics depend on the bounce being rotating (Kerr-like) vs. non-rotating (Schwarzschild-like)
- The exact relationship between rotation, chirality, and the preferred-axis signature observed in galaxy rotations

## Cross-references

- `01-cartasis.md`: the membrane where chirality dynamics happen
- `02-baby-universes.md`: what the baby inherits
- `04-supraverse.md`: lineage structure with alternating chirality
- `05-dark-sector.md`: dark matter as filter-rejected content
- `07-observational-tests.md`: specific predictions