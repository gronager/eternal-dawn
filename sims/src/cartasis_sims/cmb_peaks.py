r"""Tier 2: CMB acoustic-peak HEIGHTS from the bounce's scale-invariant spectrum.

Positions were the easy half (acoustic.py): they are set by the sound horizon, shared
with LambdaCDM. The HEIGHTS are the discriminating half -- they encode the baryon and
dark-matter densities through recombination physics, again shared with LambdaCDM, so
SCT inherits them once the bounce supplies the (near scale-invariant) primordial
input. The single most physical, falsifiable feature is the ODD/EVEN peak
alternation, driven by the baryon loading

    R = (3/4)(rho_b/rho_gamma) = (3/4)(omega_b/omega_gamma)/(1+z_rec),

which at recombination is R ~ 0.6. Baryons add inertia to the photon-baryon fluid, so
COMPRESSION peaks (1st, 3rd, 5th) are enhanced and RAREFACTION peaks (2nd, 4th)
suppressed: the effective-temperature monopole at last scattering is

    Theta_eff(k) = (1+R) Psi cos(k r_s) - R Psi,

whose extrema have magnitude (1+2R) at compressions and 1 at rarefactions. A Doppler
(velocity) term fills the troughs, Silk diffusion damps the tail, and radiation
driving boosts the first peaks. Here we build that standard semi-analytic
tight-coupling spectrum -- R and the acoustic scale are COMPUTED; the overall
amplitude, the damping scale, and a couple of O(1) shape factors are calibrated --
and show it reproduces the Planck peak heights. A Boltzmann code (CAMB/CLASS) refines
it, but the heights follow from omega_b (which SCT must deliver from the inherited
eta, Ch.4) and omega_c (the dark sector, Ch.6).
"""

from __future__ import annotations

import numpy as np

from . import acoustic as ac

# ---- computed inputs -------------------------------------------------------
OMEGA_GAMMA_HSQ = 2.47e-5      # photon density
Z_REC = ac.Z_REC


def baryon_loading(omega_b_hsq: float = 0.02237, z_rec: float = Z_REC) -> float:
    """R = (3/4)(omega_b/omega_gamma)/(1+z_rec): baryon-to-photon inertia ratio
    at recombination (~0.6). Larger R -> stronger odd/even peak alternation."""
    return 0.75 * (omega_b_hsq / OMEGA_GAMMA_HSQ) / (1.0 + z_rec)


def compression_rarefaction_ratio(omega_b_hsq: float = 0.02237) -> float:
    """First-to-second (compression/rarefaction) peak HEIGHT ratio in D_l. The crude
    tight-coupling monopole gives (1+2R)^2, but the WKB amplitude ~(1+R)^{-1/4} and
    radiation driving soften it; we use the calibrated baryon shift A_b ∝ R, ratio
    = ((1+A_b)/(1-A_b))^2, which gives ~2.2 at the Planck omega_b (vs ~5 for the
    crude monopole)."""
    A_b = baryon_shift(omega_b_hsq)
    return ((1.0 + A_b) / (1.0 - A_b)) ** 2


def baryon_shift(omega_b_hsq: float = 0.02237) -> float:
    """Effective baryon shift A_b of the acoustic zero-point, proportional to the
    loading R and calibrated so the first-to-second peak ratio is ~2.2 at the Planck
    omega_b. The monopole is cos(k r_s) - A_b: compressions ~(1+A_b), rarefactions
    ~(1-A_b)."""
    return _K_BARYON * baryon_loading(omega_b_hsq)


# ---- semi-analytic tight-coupling spectrum --------------------------------
# Computed: acoustic scale l_A and R. Calibrated O(1): the baryon-shift coefficient
# (so the alternation matches), amplitude, Silk scale, Doppler weight, driving.
_PHASE = 0.25                  # acoustic phase (l_1 ~ 0.75 l_A)
_K_BARYON = 0.305              # A_b = K_baryon * R  -> A_b ~ 0.19 at Planck omega_b
_DOPP = 0.60                   # Doppler (velocity) relative weight (fills troughs)
_ELL_D = 830.0                 # Silk damping multipole
_DAMP_P = 1.5                  # Silk damping power
_ELL_DRIVE = 200.0             # radiation-driving turn-on
_DRIVE = 2.3                   # driving boost of acoustic over Sachs-Wolfe
_SW = 0.55                     # Sachs-Wolfe plateau amplitude (rel.)
_AMP = 1.9e3                   # overall normalization [microK^2]


def _acoustic_scale_l():
    _, lA, _, _ = ac.acoustic_scale()
    return lA


def cl_tt(ell, omega_b_hsq: float = 0.02237, n_s: float = 0.965, **kw):
    """Semi-analytic CMB TT spectrum D_l = l(l+1)C_l/2pi [microK^2].

    Baryon-loaded acoustic monopole + Doppler, radiation-driving boost, Silk
    damping, and a Sachs-Wolfe plateau, on a near scale-invariant primordial tilt."""
    ell = np.asarray(ell, dtype=float)
    lA = kw.get("lA", _acoustic_scale_l())
    A_b = baryon_shift(omega_b_hsq)
    R = baryon_loading(omega_b_hsq)
    arg = np.pi * (ell / lA + _PHASE)
    monopole = np.cos(arg) - A_b                     # baryon-shifted: enhances odd peaks
    doppler = np.sin(arg) / np.sqrt(1.0 + R)         # fills the troughs
    osc = monopole**2 + _DOPP * doppler**2
    driving = 1.0 + _DRIVE / (1.0 + (_ELL_DRIVE / np.maximum(ell, 1.0)) ** 2)
    silk = np.exp(-(ell / _ELL_D) ** _DAMP_P)        # diffusion damping
    primordial = (ell / 220.0) ** (n_s - 1.0)
    sw = _SW / (1.0 + (ell / 25.0) ** 2) ** 0.6      # large-scale plateau
    return _AMP * primordial * (sw + osc * driving * silk)


def peak_heights(omega_b_hsq: float = 0.02237, **kw):
    """Return (peak_ells, peak_Dl) of the first ~5 acoustic peaks of the model."""
    ell = np.arange(2, 2400.0)
    D = cl_tt(ell, omega_b_hsq, **kw)
    # local maxima above the Sachs-Wolfe rise
    mx = (D[1:-1] > D[:-2]) & (D[1:-1] > D[2:]) & (ell[1:-1] > 150)
    pe, pd = ell[1:-1][mx], D[1:-1][mx]
    return pe[:5], pd[:5]
