r"""The ab-initio spectrum: the fermion mass matrix and W/Z/H from G, hbar, c ALONE.

This is the strict Weltformel computation (Heisenberg's 1950s-60s nonlinear-spinor program):
ONE nonlinear Dirac equation -- the Einstein-Cartan four-fermion (Hehl-Datta) soliton --
with NO inserted Yukawas, NO electroweak v, NO f_pi, NO fitted well. The only dimensionful
inputs are the three constants of Nature:

    G (Newton), hbar (Planck), c (light).

Because the four-fermion coupling IS gravity (kappa = 8 pi G / c^4), it sets its own scale:
the soliton's natural mass is the PLANCK MASS, M_Pl = sqrt(hbar c / G) ~ 1.22e19 GeV. So the
absolute numbers come out ~10^17-10^22 too large -- 'off by a lot, but ab initio', as
expected. What is genuine ab-initio CONTENT is the dimensionless STRUCTURE:

  * generations  -> the radial eigenlevels & their wavefunction overlaps (self_consistent.py,
                    NO fitted shape -- the well is sourced by the levels it binds);
  * towers       -> the charge/colour 'handles' (ab_initio_charges.py, NO fit), which give
                    quarks > leptons and force the NEUTRINO to zero;
  * bosons       -> the composite scales of the same soliton (the condensate peak ~ f_pi,
                    the vector excitation gap ~ M_V) set W, Z, and the composite Higgs.

The output is a 3x4 fermion matrix and the three boson masses, all in GeV, all Planck-scaled.
The gap to the observed values is the SCALE (Planck-vs-electroweak) hierarchy -- the density
question (Part II target L5): whether the particle-scale binding density differs from the
cosmological rho_C, bringing M_Pl down to the electroweak scale. That single ratio is the
one number the laptop does NOT produce; everything else is turned out by the crank. Honest
1960s-style: turn the handle, read off the matrix, see the structure, own the scale gap.
"""

from __future__ import annotations

import numpy as np

from . import self_consistent as sc
from .self_consistent import _solve_levels
from . import ab_initio_charges as ai

# --- the only dimensionful inputs: G, hbar, c (SI) ---
G_NEWTON = 6.67430e-11        # m^3 kg^-1 s^-2
HBAR = 1.054571817e-34        # J s
C_LIGHT = 2.99792458e8        # m s^-1
GEV_PER_JOULE = 1.0 / 1.602176634e-10

# towers and generation labels (3 generations x 4 towers)
TOWERS = ["up-quark", "down-quark", "charged-lepton", "neutrino"]
GEN_LABELS = {
    "up-quark": ["u", "c", "t"],
    "down-quark": ["d", "s", "b"],
    "charged-lepton": ["e", "mu", "tau"],
    "neutrino": ["nu1", "nu2", "nu3"],
}
OBSERVED_GEV = {  # PDG, GeV
    "up-quark": [2.16e-3, 1.273, 172.57],
    "down-quark": [4.7e-3, 93.5e-3, 4.183],
    "charged-lepton": [0.511e-3, 105.66e-3, 1.77693],
    "neutrino": [1e-11, 9e-11, 5e-10],   # ~ order 0.01-0.05 eV, illustrative
}


def planck_mass_GeV(G=G_NEWTON, hbar=HBAR, c=C_LIGHT):
    """The natural scale of the Einstein-Cartan four-fermion soliton: M_Pl c^2 = sqrt(hbar
    c^5 / G), in GeV. This is the ONE scale, and it comes from G, hbar, c with nothing else."""
    E_planck_joule = np.sqrt(hbar * c**5 / G)
    return E_planck_joule * GEV_PER_JOULE


def soliton_structure(g=4.0, m_sigma=1.0, n_gen=3):
    """The dimensionless soliton, NO fitted shape: solve the self-consistent four-fermion
    well and read off (a) the level energies E_n (generations) and (b) the wavefunction
    overlaps with the condensate O_n (the configurational coupling). Returns dimensionless
    arrays -- the structure that modulates the Planck scale."""
    out = sc.solve_soliton(m0=1.0, g=g, m_sigma=m_sigma, n_fermions=1,
                           n_levels=max(4, n_gen + 1))
    r = out["r"]
    h = r[1] - r[0]
    _, u = _solve_levels(out["M"], r, h, n_gen)
    sig = out["sigma"]
    overlaps = []
    for n in range(n_gen):
        un = u[:, n]
        overlaps.append(float(np.trapezoid(un**2 * sig, r) / np.trapezoid(un**2, r)))
    return {
        "levels": np.asarray(out["E"][:n_gen]),
        "overlaps": np.asarray(overlaps),
        "soliton_mass": float(out["mass"]),
        "f_pi": float(out["f_pi"]),       # condensate peak (dimensionless)
        "M_V": float(out["M_V"]),         # vector-excitation gap (dimensionless)
        "converged": bool(out["converged"]),
    }


def fermion_matrix_GeV(g=4.0, m_sigma=1.0):
    """The 3x4 fermion mass matrix, ab initio (G, hbar, c + charges only). Each entry is

        m(T, n) = M_Pl * handle_norm(T) * overlap_norm(n) ,

    with the Planck mass as the scale, the charge/colour handle as the tower factor (no fit,
    zero for the neutrino), and the soliton overlap as the generation factor (no fit). All
    entries are ~Planck-scaled -- off from observation by the hierarchy, by construction."""
    M_Pl = planck_mass_GeV()
    st = soliton_structure(g, m_sigma)
    ov = st["overlaps"]
    # the heaviest fermion is the most-bound (largest-overlap) state; assign generations by
    # overlap magnitude (gen III = largest overlap, gen I = smallest). This is the physical
    # identification, NOT a fit -- the ground state is the heaviest rung.
    ov_norm = np.sort(ov) / np.sort(ov).max()                # ascending: [genI .. genIII], max=1
    handles = {t: ai.condensate_handle(t) for t in TOWERS}
    hmax = max(handles.values())
    matrix = {}
    for t in TOWERS:
        hn = handles[t] / hmax                                # tower factor, max=1
        matrix[t] = M_Pl * hn * ov_norm
    return {"matrix_GeV": matrix, "M_Pl_GeV": M_Pl,
            "overlaps": ov, "gen_factor": ov_norm, "handles": handles, "structure": st}


def boson_masses_GeV(g=4.0, m_sigma=1.0):
    """W, Z, Higgs from the SAME soliton's composite scales, ab initio. The condensate peak
    plays f_pi and the vector gap plays M_V; with the (input, not fitted) electroweak mixing
    these give the boson scale = M_Pl * (dimensionless composite). Planck-scaled like the
    fermions. The Higgs is the composite scalar (the sigma itself)."""
    M_Pl = planck_mass_GeV()
    st = soliton_structure(g, m_sigma)
    fpi, MV = st["f_pi"], st["M_V"]
    sin2 = 0.23121
    # m_W ~ (1/2) g_weak F with F ~ fpi*M_Pl; here we expose the composite scale F and the
    # custodial structure, all in Planck units (g_weak is a gauge coupling, not a fit).
    F = fpi * M_Pl                                            # composite EW scale (GeV)
    g_weak = np.sqrt(4 * np.pi / 128.0) / np.sqrt(sin2)
    mW = 0.5 * g_weak * F
    mZ = mW / np.sqrt(1 - sin2)
    mH = MV * M_Pl                                            # composite scalar ~ vector gap
    return {"W": mW, "Z": mZ, "H": mH, "composite_F_GeV": F, "M_Pl_GeV": M_Pl}


def hierarchy_gap():
    """The one number the laptop does NOT produce: the ratio of the ab-initio (Planck) scale
    to the observed electroweak scale -- the gauge hierarchy = the density question (L5).
    Returns M_Pl / (top mass), the factor by which the whole matrix is too heavy."""
    M_Pl = planck_mass_GeV()
    return M_Pl / OBSERVED_GEV["up-quark"][-1]
