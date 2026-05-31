# What sets the Cartasis density — and does chirality add a counter-pressure?

Working note, prompted by: is 10^50 kg/m^3 fixed, or could other limits change
the number and the dynamics — e.g. a chirality counter-pressure unknown in
Cartan's day? Short answer: the number is a composition-dependent *estimate*, not
a postulate, and yes — there is a plausible additional, parity-violating chiral
counter-pressure that (a) is genuinely post-Cartan physics and (b) is exactly the
P-odd ingredient Chapter 4 said a real filter would need.

## 1. rho_C is an estimate with moving parts

The bounce threshold is where torsion repulsion (~ spin density squared) overtakes
gravity (~ energy density):

    rho_C ~ c^4 / (G hbar^2 (s/rho)^2),    with s/rho ~ hbar/m_f
        =>  rho_C ~ m_f^2 c^4 / (G hbar^4).

The 10^50 figure used the **nucleon** mass for m_f ("one spin per nucleon"). But
at 10^50 kg/m^3 matter is far beyond quark deconfinement — there are no nucleons.
The spin is carried by quarks and leptons, whose relevant (current) masses are
much smaller, so s/rho is *larger* and:

> **rho_C scales as m_f^2.** If light fermions carry the spin, the bounce can
> trigger at a substantially *lower* density than 10^50 — possibly several orders
> lower. The number is composition-dependent (open question Q1/Q2), and 10^50 is
> best read as "nucleonic-spin" upper estimate.

Plus the usual coefficient ambiguity (factors 10–100 from the precise
Einstein–Cartan coefficients).

## 2. The pressures on the way in (the equation of state)

Three things resist or shape the collapse before/at the bounce:

- **Degeneracy (Pauli) pressure.** Holds up white dwarfs and neutron stars, but
  ultra-relativistic degeneracy pressure scales like the energy density's own
  contribution to gravity (rho^{4/3}) — it sets the inbound equation of state
  (our Tier-1 sweep w, tending stiff), but it does *not* halt collapse. That is
  why neutron stars above the TOV limit become black holes rather than stalling.
- **QCD phase structure.** Through the relevant range the matter passes
  deconfinement and likely color-superconducting (quark-pairing) phases. Pairing
  is *attractive* and could soften the EoS and partially cancel spin density
  (paired spins) — nudging rho_C and the bounce sharpness. Unmodelled so far.
- **Torsion (the axial-current-squared term).** This is the actual repulsion that
  wins because it grows as rho^2, faster than gravity's rho. It *is* a spin /
  chirality counter-pressure — see below.

So degeneracy sets w going in; torsion provides the bounce; QCD phases modulate
both. None of these is in the bare 10^50 estimate.

## 3. The chirality counter-pressure — and why it is "post-Cartan"

Cartan formulated torsion in 1922–25. The **axial/chiral structure of matter**
— parity violation (1956), the V–A weak current, the axial current
j5 = psi-bar gamma^5 gamma^mu psi that torsion actually couples to — is all
post-1950s. So the bounce repulsion is, at root, a *chiral* effect that Cartan
could not have known he was writing down. The user's instinct is well-placed.

Beyond the minimal Hehl–Datta (axial)^2 term, at rho_C the matter is a hot,
dense, relativistic **chiral plasma**, where modern anomaly physics applies and
adds genuinely chirality-dependent pressures:

- **Chiral imbalance free energy.** A net chirality (n_L != n_R, i.e. a chiral
  chemical potential mu5) costs free energy; restoring balance releases it. This
  is a real chirality-dependent pressure/relaxation term.
- **Chiral magnetic effect (CME):** a current j ∝ mu5 B along a magnetic field;
- **Chiral vortical effect (CVE):** a current along the **vorticity** (rotation).

CME and CVE are **parity-odd** and tie chirality to magnetic fields and to
*rotation*. The anomaly can even change n_L - n_R (via gauge-field topology /
sphalerons), so chirality is not simply conserved at these temperatures.

## 4. The payoff: this is where a real filter (and maybe the asymmetry) lives

Chapter 4 showed the minimal Einstein–Cartan term is C/P/T-even, so it cannot
flip or filter. We said a genuine filter needs a parity-violating ingredient and
pointed at exotic gravity (Holst/Nieh–Yan). But **chiral-plasma anomaly physics
is parity-violating Standard-Model physics** — it may be the P-odd ingredient we
needed, with no new gravity required:

- The **CVE couples chirality to vorticity** — i.e. to the parent black hole's
  **rotation/spin**. A rotating (Kerr) bounce therefore has a built-in axis along
  which chiral transport is biased: a *derived* chirality filter, tied directly
  to the spin axis of the axis-of-evil / galaxy-spin story.
- Combined with the existing Standard-Model CP violation and the
  out-of-equilibrium bounce, this is the Sakharov-condition machinery operating
  in the bounce background — a route to a *real* asymmetry (a source, not just
  the mean-field amplifier of Chapter 4), without postulating the CPT mirror.

So the "chirality counter-pressure" question reconnects three loose ends at once:
it can move rho_C (composition), it modifies the bounce dynamics (chiral EoS),
and — most importantly — it is a concrete, modern, computable candidate for the
filter/baryogenesis gap, naturally aligned with the parent's spin axis.

## 5. To compute

1. **rho_C(composition):** recompute the threshold with quark/lepton spin
   carriers and the QCD-phase EoS; bracket how far below 10^50 it can sit.
2. **Chiral EoS through the bounce:** add chiral imbalance pressure to the Tier-1
   integrator; see how it shifts a_min, the sharpness, and rho_C.
3. **CVE filter in the rotating bounce:** estimate the handedness bias from
   chiral-vortical transport along the spin axis at rho_C; compare to the
   dark-matter ratio f ~ 1/6 and the baryon asymmetry eta ~ 6e-10.
4. Feed the verdict back into Chapter 4 (filter = anomalous chiral transport,
   P-odd, spin-axis-aligned) and the axis story (Chapter 8).
