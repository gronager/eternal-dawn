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

**What it establishes (robust):**
- The torsiton **binds** — a clean effective-mass plateau at m_N > m_π > 0, non-perturbatively, with
  no mean-field approximation. (A wrong-sign contraction would give a negative/junk correlator; this
  is clean and correctly ordered.)
- **Constituent counting.** At heavy quark mass m_N/m_π → 3/2 (a baryon is three quarks, a pion two);
  we sit just above it (1.635 at the heaviest point). This is the lattice confirming the mean-field
  result **M_torsiton ≈ N_c·M** (Appendix C, `kahana_ripka`) — the two halves of the program agree.
- **Chiral trend.** m_π² is linear in the bare mass (GMOR), giving a critical mass m_crit ≈ −1.12;
  as the quark lightens, m_N/m_π **rises** (1.635 → 1.871, the increments accelerating) — the mass
  going condensate-dominated. Rough chiral extrapolation: m_N(chiral) a ≈ 1.4, i.e. m_N/√σ ≈ 2.8
  (β=5.6 → √σ a ≈ 0.5).

**What it does NOT establish (the honest caveats):** this is **quenched** (no sea quarks), **coarse**
(β=5.6), and **small** (12³). All three bias the number — the chiral m_N/√σ ≈ 2.8 here vs a real
baryon's ≈ 2.1. The pilot proves the *object and its qualitative physics*, not the physical mass. The
quenched scan walls near m_crit where exceptional configurations (near-zero Wilson modes, unsuppressed
without a sea) make the plateau wild.

## L4 production — dynamical (in progress)

The physical mass, the true (fluctuating) condensate, and any excited rungs (candidate further
generations) need the **fermion determinant**. `generate_dynamical` (SU(3) N_f=2 fundamental Wilson
HMC) + `run/07_dynamical_torsiton.sh` set this up; the same `measure_baryon` spectroscopy runs on the
dynamical configs. Leadership-class cost (the CG sea force is dear per trajectory) — start small,
watch acceptance and dH, scale up.

## L1 — confinement (status)

The pure-gauge area law / string tension pipeline (`measure_potential`, `run/01`) is in hand and
audited; it gives the scale √σ that turns the torsiton mass into a cutoff-independent ratio. The
remaining L1 piece (a dual order parameter for the confinement *mechanism*) is open.
