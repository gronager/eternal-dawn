r"""Tier 2, toward the acoustic peaks: sound horizon, peak positions, schematic C_l.

The bounce supplies a near-scale-invariant PRIMORDIAL spectrum (primordial.py).
The CMB acoustic peaks are then a RECOMBINATION-era feature -- baryon--photon
sound waves frozen at last scattering -- and Supraverse Cartasis Theory shares the
entire post-bounce hot Big Bang with LambdaCDM, so it inherits the peaks: same
recombination physics acting on the same plasma. (The "CMB = parent Hawking
radiation" identification of Ch.6 is about the ORIGIN of that thermal plasma at the
bounce/horizon surface, not a replacement for recombination; the peaks come from
standard pre-recombination acoustics either way.)

The acoustic scale is set by the comoving sound horizon at recombination

    r_s = Integral_{z_rec}^{inf} c_s(z)/H(z) dz,   c_s = c / sqrt(3(1+R)),
    R(z) = (3 rho_b)/(4 rho_gamma) = (3 Omega_b)/(4 Omega_gamma) * 1/(1+z),

and the comoving distance to last scattering D_C = Integral_0^{z_rec} c/H dz, giving
the acoustic angular scale theta_s = r_s/D_C and peak multipoles l_n ~ n pi/theta_s
(the first peak shifted to l_1 ~ 220 by the well-known ~0.75 phase factor). With
Planck cosmology this returns r_s ~ 145 Mpc and l_1 ~ 220 -- the observed value.
A schematic C_l (scale-invariant primordial x a damped acoustic transfer) shows the
peak series; it is illustrative, not a Boltzmann computation.
"""

from __future__ import annotations

import numpy as np

from . import constants as k

Z_REC = 1089.8                 # redshift of recombination (Planck)
T_CMB = 2.7255                 # K
# photon density parameter Omega_gamma h^2 from T_cmb
OMEGA_GAMMA_HSQ = 2.47e-5      # ~ (T/2.725)^4 * 2.47e-5
N_EFF_NU = 3.046


def _hubble(z, Om, OL, h, Or):
    """H(z) [s^-1] for a flat LambdaCDM background (matter+radiation+Lambda)."""
    H0 = h * 1.0e5 / k.Mpc      # 100h km/s/Mpc in s^-1
    return H0 * np.sqrt(Om * (1 + z) ** 3 + Or * (1 + z) ** 4 + OL)


def _R_baryon(z, Ob_hsq):
    """Baryon-to-photon momentum density ratio R = 3 rho_b/(4 rho_gamma)."""
    Og_hsq = OMEGA_GAMMA_HSQ * (1.0 + 0.2271 * N_EFF_NU)   # photons (+nu in c_s? no)
    return (3.0 * Ob_hsq) / (4.0 * OMEGA_GAMMA_HSQ) / (1.0 + z)


def sound_horizon(Om=0.3153, OL=0.6847, h=0.6736, Ob_hsq=0.02237,
                  z_rec=Z_REC, n=20000):
    """Comoving sound horizon r_s at recombination [Mpc]."""
    Or = OMEGA_GAMMA_HSQ * (1.0 + 0.2271 * N_EFF_NU) / h**2
    z = np.linspace(z_rec, z_rec + 1.0e5, n)
    cs = k.c / np.sqrt(3.0 * (1.0 + _R_baryon(z, Ob_hsq)))   # m/s
    integ = cs / _hubble(z, Om, OL, h, Or)                   # m
    return float(np.trapezoid(integ, z) / k.Mpc)             # Mpc


def comoving_distance_lss(Om=0.3153, OL=0.6847, h=0.6736, z_rec=Z_REC, n=20000):
    """Comoving distance to last scattering D_C [Mpc]."""
    Or = OMEGA_GAMMA_HSQ * (1.0 + 0.2271 * N_EFF_NU) / h**2
    z = np.linspace(0.0, z_rec, n)
    integ = k.c / _hubble(z, Om, OL, h, Or)
    return float(np.trapezoid(integ, z) / k.Mpc)


def acoustic_scale(**kw):
    """Acoustic angular scale theta_s = r_s/D_C (radians) and l_A = pi/theta_s."""
    rs = sound_horizon(**{q: kw[q] for q in kw if q in
                          ("Om", "OL", "h", "Ob_hsq", "z_rec")})
    dc = comoving_distance_lss(**{q: kw[q] for q in kw if q in
                                  ("Om", "OL", "h", "z_rec")})
    theta = rs / dc
    return theta, np.pi / theta, rs, dc


def peak_multipoles(n_peaks=5, phase=0.75, **kw):
    """Peak positions l_n ~ (n - 1 + phase) pi/theta_s; phase~0.75 puts l_1~220."""
    _, lA, _, _ = acoustic_scale(**kw)
    nn = np.arange(1, n_peaks + 1)
    return (nn - 1.0 + phase) * lA


def schematic_cl(ell, n_s=0.965, l_damp=1400.0, phase=0.75, **kw):
    """Schematic CMB TT power D_l = l(l+1)C_l/2pi (arbitrary norm): an oscillatory
    acoustic transfer (peaks every l_A, first at ~0.75 l_A) x Silk damping x a
    large-scale turnover, with the scale-invariant primordial input nearly flat.
    Illustrative only -- not a Boltzmann-code spectrum. Peaks sit at peak_multipoles."""
    ell = np.asarray(ell, dtype=float)
    _, lA, _, _ = acoustic_scale(**kw)
    primordial = ell ** (n_s - 1.0)                    # ~ scale-invariant input
    osc = (0.55 + 0.45 * np.cos(2.0 * np.pi * (ell / lA - phase))) ** 2
    damp = np.exp(-(ell / l_damp) ** 1.8)              # Silk diffusion damping
    rise = (ell / 30.0) ** 2 / (1.0 + (ell / 30.0) ** 2)   # large-scale turnover
    return primordial * osc * damp * rise
