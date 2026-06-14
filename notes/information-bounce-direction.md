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
