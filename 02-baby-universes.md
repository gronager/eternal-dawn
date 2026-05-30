# 02 — Baby Universes

## What forms after the bounce

When matter reaches Cartasis density and the torsion repulsion overwhelms gravity, the bouncing matter doesn’t return to the parent’s exterior (the parent’s horizon is now in its past, geometrically). Instead, the bounce dynamics generate a new spacetime region — a baby universe — with its own time direction and spatial extent.

From inside, the baby universe looks like a standard hot-Big-Bang cosmology:

- Begins at maximum density (the bounce moment)
- Expands rapidly (torsion-driven, qualitatively similar to inflation)
- Cools as it expands
- Transitions to standard radiation-dominated and matter-dominated phases
- Eventually shows large-scale structure formation

From outside (parent’s perspective), the baby is not visible as a spatial structure. The parent sees a black hole at the location where the collapse happened. The baby exists “in the future of the bounce” — its spatial extent is in directions orthogonal to the parent’s interior geometry.

## The Tardis effect

A consequence of the geometry: the parent black hole, viewed from outside, has a specific Schwarzschild radius. The baby universe inside has its own cosmological extent, which can be vastly larger than the parent’s horizon would suggest.

For our universe: total mass-energy is ~$10^{53}$ kg. At Cartasis density, this fits in ~$10^3$ m³ — less than half an Olympic swimming pool. But from inside, we observe a universe ~46 billion light-years in radius (observable horizon).

This is not a contradiction. The two volumes are measured in different spacetime regions. The “small volume” measures the matter at the bounce moment in parent coordinates. The “large volume” measures the cosmological extent in baby coordinates after expansion. They’re not the same quantity.

This is also not unique to Einstein-Cartan. Even in standard GR, the spatial volume inside a black hole’s horizon grows with time (Christodoulou-Rovelli) and becomes essentially unbounded. The bounce framework reinterprets this growth as cosmological expansion of the baby rather than approach to a singularity.

## Inflation substitute

The early baby universe undergoes rapid expansion because torsion repulsion is at its maximum just after the bounce, when density is still very high. As the baby expands and density drops, the torsion repulsion fades. Gravity reasserts dominance, expansion decelerates.

This phase has the right qualitative features to substitute for the inflationary epoch of standard cosmology:

- Brief, very rapid expansion
- Ends naturally when conditions change (density drops below Cartasis threshold)
- No separate inflaton field required
- Sets up homogeneous and approximately isotropic conditions for subsequent cosmological evolution

Quantitative work needed: derive the power spectrum of perturbations produced by torsion-driven expansion and compare to inflation’s predictions and to CMB observations. This is the central observational test of the bounce framework (see `07-observational-tests.md` and `08-simulation-plan.md`).

## What the baby inherits from the parent

Several quantities are transmitted across Cartasis:

**Energy-momentum.** Conservation laws ensure the baby’s total energy-momentum equals what fell into the bounce. This sets the baby’s mass-energy budget.

**Spin and angular momentum.** The parent’s net angular momentum (if rotating) is transmitted through the bounce. A rotating parent (Kerr-like) produces a baby with a preferred axis — possibly the source of cosmological anisotropies like the galaxy rotation asymmetry observed in JWST data (Shamir 2025).

**Chirality bias.** The parent’s net matter-antimatter asymmetry is transmitted, possibly with flip and/or filter modifications (see `03-chirality-flip-filter.md`).

**Quantum entanglement.** Particles in the parent’s near-horizon region remain entangled with their Hawking radiation partners. This entanglement persists across the bounce, connecting the baby’s interior with the parent’s exterior — see `05-dark-sector.md`.

**Conserved charges.** Baryon number, lepton number, electric charge (if conserved at trans-EW temperatures) flow through.

What is **not** transmitted:

- Structure (atoms, molecules, nuclei) — all dismantled at Cartasis density
- Specific particle identities — only conserved currents survive
- Information about the parent’s “before” — the bounce is the baby’s past boundary

## The persistent channel

As discussed in `01-cartasis.md`, the bounce is not a single event but an ongoing interface. While the parent black hole continues to accrete matter, that matter falls in, reaches Cartasis, and crosses into the baby. From the baby’s perspective, this means matter continues to be injected throughout the parent’s lifetime.

The injection geometry is non-trivial: new matter doesn’t appear at a specific spatial location in the baby’s 3D space. It appears “at the bounce membrane,” which from the baby’s perspective is the cosmic horizon at maximum lookback. So injected matter manifests as gravitational influence from the cosmic boundary, not as visible particles appearing in our 3D space.

This is the basis for the dark matter interpretation in `05-dark-sector.md`.

## Lifecycle

A baby universe’s lifetime is bounded by its parent black hole’s Hawking evaporation timescale. For a parent of mass $M$:

$$\tau_H \sim \frac{5120 \pi G^2 M^3}{\hbar c^4}$$

For a stellar-mass parent: $\tau_H \sim 10^{67}$ years.
For a supermassive parent ($10^9 M_\odot$): $\tau_H \sim 10^{94}$ years.
For our universe’s mass ($\sim 10^{53}$ kg): $\tau_H \sim 10^{145}$ years.

These vastly exceed any internal cosmological timescale of the baby (~$10^{10}$ years for matter-radiation transitions, ~$10^{100}$ years for stellar evaporation). So the baby’s internal physics plays out fully before parent-evaporation becomes relevant.

But the baby’s lifetime is finite. At some point, the parent’s mass is mostly radiated away, the bounce membrane shrinks, and the baby’s connection to the parent fades. What happens to the baby at that point is unclear — possibly internal physics continues, possibly the baby’s spacetime structure also dissolves.

## Generation of black holes within the baby

The baby’s internal physics, if similar to ours, produces stellar collapse, supernovae, and astrophysical black holes. Each such black hole, in the bounce framework, is itself a potential progenitor of further baby universes — grandchildren in the supraverse lineage.

Stellar mass black holes (~10 solar masses) have Cartasis-density throat regions that produce relatively small grandchildren. Supermassive black holes (~$10^9 M_\odot$) produce larger grandchildren, potentially comparable to our own universe in scale.

This cascading structure means each viable universe produces many descendants, and the supraverse fills with universe-trees rooted in OG bounces (see `04-supraverse.md`).

## Cross-references

- `01-cartasis.md`: the membrane that the baby emerges from
- `03-chirality-flip-filter.md`: matter-antimatter dynamics across the bounce
- `04-supraverse.md`: global structure of nested baby universes
- `05-dark-sector.md`: observational consequences of the persistent channel
- `06-time-arrow.md`: the bounce as the baby’s low-entropy boundary
- `08-simulation-plan.md`: numerical evolution of bounces