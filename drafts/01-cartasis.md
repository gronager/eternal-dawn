# 01 — Cartasis

## The bounce membrane

**Cartasis** is the 3D hypersurface in spacetime where torsion-mediated repulsion equals gravitational attraction, and across which the geometry transitions from collapsing-parent to expanding-baby. The name captures both the surface-as-membrane and the moment-as-genesis aspects: it is simultaneously a spatial boundary in supraverse coordinates and a temporal beginning in baby-universe coordinates.

## Physics of the bounce

### The competing terms

In Einstein-Cartan with fermionic matter, the effective stress-energy at high density has two relevant contributions:

**Gravitational term** (attractive at all densities):
$$T_{\text{grav}} \sim \rho c^2$$
scaling linearly with mass-energy density.

**Torsion term** (repulsive at high spin density):
$$T_{\text{tors}} \sim G \hbar^2 s^2 / c^4$$
scaling as the square of spin density $s$.

At ordinary densities, $T_{\text{tors}} / T_{\text{grav}} \sim 10^{-40}$, negligible. As matter compresses, both terms grow, but spin density grows faster (because compression packs more fermions into each volume, each contributing $\hbar/2$). The quadratic scaling means $T_{\text{tors}}$ overtakes $T_{\text{grav}}$ at sufficient compression.

### The critical density

Setting $T_{\text{grav}} = T_{\text{tors}}$ and solving for density:

$$\rho_C \sim \frac{c^4}{G \hbar^2 (s/\rho)^2}$$

For typical fermionic matter where $s/\rho \sim \hbar/m_n$ (one fermion-spin per nucleon-mass), this gives:

$$\rho_C \sim \frac{m_n^2 c^4}{G \hbar^4} \sim 10^{50} \text{ kg/m}^3$$

This is 33 orders of magnitude beyond nuclear density and 46 orders of magnitude below Planck density. The exact value depends on the fermion content (different masses give different ratios) and on detailed coefficients in the Einstein-Cartan Lagrangian.

### What’s actually at Cartasis density

At $\rho_C$, conventional particle descriptions break down. The matter is compressed beyond quark deconfinement, beyond electroweak symmetry restoration, beyond any regime where the Standard Model’s standard quasiparticle excitations are well-defined. What exists is **continuous field content**: energy-momentum density, spin density, chiral current density. Individual particle identities don’t survive.

What crosses Cartasis is therefore not “particles” but **conserved currents**: total energy-momentum, total spin, total chiral charge, and (subject to electroweak anomaly considerations) baryon and lepton numbers.

This has implications for chirality dynamics (see `03-chirality-flip-filter.md`) and for what kind of structures can be transmitted through a bounce (none — structures must be rebuilt from elementary field content on the other side).

## Geometry of the bounce

### The radial-temporal exchange

Inside a Schwarzschild black hole’s event horizon, the radial coordinate $r$ becomes timelike. Future-directed motion is radially inward. This is standard GR.

At Cartasis density, the bounce dynamics reverse this: the torsion repulsion pushes matter outward in $r$, but since $r$ is timelike inside the horizon, “outward in $r$” doesn’t take matter back through the horizon — it takes matter into the *future* of the bounce.

The future of the bounce, in baby universe coordinates, is the cosmological future of the new spacetime region. The radial direction $r$ on the parent’s side smoothly transitions into the time direction $\tau$ on the baby’s side at the bounce hypersurface.

### Junction conditions

The metric and connection must be continuous across the bounce. Specifically:

- Induced metric on the Cartasis hypersurface matches from both sides
- Extrinsic curvature matches (modified Israel-Darmois conditions accommodating torsion)
- Stress-energy and spin current flow continuously through the surface

These junction conditions are the technical specification of what “the bounce is continuous” means. They’re a system of equations relating parent-side and baby-side field values at the membrane.

### Bounce as persistent interface, not single event

Critical point: in a black hole that continues to accrete matter from its parent universe’s environment, the bounce is **ongoing**. New matter falls in, compresses, reaches Cartasis density, and joins the bounce. The Cartasis hypersurface is therefore not a single moment but a persistent interface that exists for as long as the parent black hole exists.

From the baby universe’s interior, this means matter continues to be injected across the Cartasis membrane throughout the parent’s lifetime — see `05-dark-sector.md` for the dark matter interpretation this enables.

## The geometric identification: parent’s horizon = baby’s cosmic horizon

From the parent’s exterior, Cartasis sits inside the event horizon. From the parent’s interior, it’s the surface where infalling matter is currently reaching critical density. From the baby’s interior, looking backward in time, Cartasis is **the boundary of the past at maximum cosmic lookback**.

The Big Bang membrane (from the baby’s perspective) and the parent’s event horizon region (from the parent’s perspective) are the same surface, viewed from different sides. This is the geometric basis for the CMB-as-Hawking-radiation identification (see `05-dark-sector.md`).

## Calculation targets

To make Cartasis quantitative, the following calculations are needed:

1. **Precise $\rho_C$ from the Einstein-Cartan Lagrangian.** Existing estimates vary by factors of 10-100 depending on coefficient choices. A clean calculation would pin the value.
1. **Bounce timescale.** How long, in proper time, does the bounce process take? This sets the timescale for the baby’s “Big Bang” phase.
1. **Bounce efficiency.** What fraction of parent’s infalling mass-energy is transmitted to the baby’s expanding phase, vs. retained in some intermediate region, vs. scattered back?
1. **Chirality selectivity.** Does the bounce treat the two chiralities symmetrically (no filter), with mild preference (weak filter), or with strong preference (strong filter)? See `03-chirality-flip-filter.md`.
1. **Spin coherence requirements.** Does the bounce require correlated spins, or does randomized spin density suffice? Existing analysis suggests random spins work because the $s^2$ scaling is positive regardless of orientation, but this needs verification.

## Open questions specific to Cartasis

- **Is $\rho_C$ universal?** Or does it depend on the matter content (different fermion species, different bound state effects)?
- **What happens at the boundary of Cartasis density** — is there a sharp transition or a gradient region?
- **How does Cartasis interact with rotation?** Kerr-like bouncing black holes vs. Schwarzschild-like.
- **Does Cartasis preserve baryon number?** Or does the electroweak anomaly at trans-EW temperatures violate B explicitly during the bounce?
- **What’s the relationship between Cartasis and quantum gravity proper?** Cartasis is in a regime where Einstein-Cartan is presumably still classical, but full quantum gravity effects might modify the picture at higher densities.

## Cross-references

- `00-axioms.md`: A4 (Cartasis density exists)
- `02-baby-universes.md`: what the bounce produces
- `03-chirality-flip-filter.md`: chirality dynamics at the membrane
- `05-dark-sector.md`: observational consequences of the membrane
- `08-simulation-plan.md`: how to numerically model the bounce