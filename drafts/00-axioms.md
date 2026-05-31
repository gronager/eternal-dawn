# 00 — Axioms

## Minimal foundation

The framework rests on four claims, the first three of which are essentially established physics and the fourth of which is a minimal natural extension.

### A1. Quantum mechanics applies

The supraverse is described by quantum field theory. Vacuum states have nonzero fluctuations governed by the uncertainty principle. Over infinite spacetime, fluctuations of arbitrary magnitude occur with probability one.

### A2. Special relativity applies

Spacetime has Lorentzian signature. Causal structure is determined by light cones. CPT symmetry holds for any local Lorentz-invariant QFT.

### A3. Einstein-Cartan gravity applies

The connection on spacetime is allowed to be asymmetric. The antisymmetric part (torsion) couples to fermion spin via the axial current. In the absence of fermionic matter, the theory reduces to standard general relativity. With fermions, torsion produces an effective interaction between spins that scales as the square of the spin density, and is therefore always positive (repulsive at high spin density) regardless of relative spin orientation.

The Einstein-Cartan Lagrangian, in shorthand:

$$\mathcal{L}*{EC} = \frac{1}{2\kappa}(R - 2\Lambda) + \mathcal{L}*{\text{matter}}(\psi, \nabla\psi, g, T)$$

where $R$ is the curvature scalar built from the full (torsion-including) connection, and the covariant derivative $\nabla$ acts on fermionic fields $\psi$ with torsion contributions.

### A4. There exists a Cartasis density

At sufficiently high mass-energy density, torsion-mediated repulsion exceeds gravitational attraction. The critical density (Cartasis density) is approximately:

$$\rho_C \sim \frac{m_P^4}{\hbar^3} \cdot \alpha_{EC}$$

where $m_P$ is the Planck mass and $\alpha_{EC}$ is a dimensionless Einstein-Cartan coupling. Estimates place this at $\rho_C \sim 10^{50}$ kg/m³, far below Planck density ($10^{96}$ kg/m³) but well above nuclear density ($10^{17}$ kg/m³).

This is calculable from A1-A3 in principle; A4 is the explicit acknowledgment that this density exists and matters.

## What we do NOT assume

- No Big Bang as a singular event
- No initial conditions for the universe
- No fine-tuned cosmological constant or inflaton field
- No special low-entropy past as a separate postulate
- No additional dimensions beyond 4
- No additional forces or fields beyond Standard Model + gravity + torsion
- No anthropic principle as a separate axiom (it emerges from the framework)

## What follows automatically

From A1-A4 plus infinity (infinite spacetime volume), the following are consequences rather than assumptions:

1. **Spontaneous OG bounce nucleation.** Rare vacuum fluctuations reach Cartasis density. By A1 plus infinity, this happens somewhere with probability one. By A4, such fluctuations bounce rather than collapsing to singularities.
1. **Baby universe production.** The bounce geometry produces new spacetime regions (developed in 02).
1. **CPT-symmetric supraverse.** Vacuum fluctuations are statistically symmetric in matter/antimatter content. Specific fluctuations have specific biases, but the supraverse-wide distribution is symmetric.
1. **Thermodynamic equilibrium configuration.** Gravitational entropy considerations (high entropy = clumped/black-hole-rich; low entropy = smooth/uniform) drive the supraverse toward foam-of-universes as its equilibrium configuration.
1. **Local arrow of time.** Each universe has a low-entropy bounce in its past and a high-entropy future, producing a local thermodynamic time-arrow.
1. **Observer placement statistics.** Observers exist in viable universes near the peak of the supraverse population distribution.
1. **Information preservation.** Hawking radiation entangles parent exterior with baby interior. The full supraverse is unitary even though individual universes appear to lose information.

## On the role of infinity

The framework requires the supraverse to be effectively infinite — at minimum, large enough that rare fluctuations to Cartasis density occur. “Effectively infinite” in this context means: large enough relative to the characteristic timescales of vacuum fluctuation that the relevant fluctuations are statistically guaranteed to occur.

For a finite supraverse, some fluctuations might never occur and the framework would have to specify which ones do. This introduces an arbitrary cutoff that we want to avoid. The natural assumption is true infinity, which removes this arbitrariness at the cost of being philosophically heavier.

The supraverse being infinite does not require it to be temporally infinite in some absolute sense — it requires only that enough “supraverse time” exist for the relevant condensation to occur. Whether the supraverse has a beginning at the supraverse level is a question we don’t currently need to answer; we only need it to be old enough that condensation has reached equilibrium.

## On the role of Einstein-Cartan

Einstein-Cartan is the minimal natural extension of GR. It’s what you get by dropping the artificial restriction that the connection must be symmetric. The connection’s antisymmetric part (torsion) is mathematically natural and couples to a physical quantity that exists in nature (fermion spin). At ordinary densities, torsion effects are utterly negligible (~10⁻⁴⁰ corrections to GR). At Cartasis density, they dominate.

We don’t postulate Einstein-Cartan as a new theory; we recognize it as the geometric framework that should have been used all along, with standard GR as the low-density limit. This is not an additional assumption beyond standard physics — it’s the removal of an assumption (symmetry of the connection).

## Cross-references

- `01-cartasis.md`: detailed treatment of A4 and the bounce mechanism
- `02-baby-universes.md`: derivation of consequence 2 (baby universe production)
- `04-supraverse.md`: derivation of consequences 3 and 4 (foam structure, equilibrium)
- `06-time-arrow.md`: derivation of consequence 5 (arrow of time)
- `07-observational-tests.md`: how to check these predictions