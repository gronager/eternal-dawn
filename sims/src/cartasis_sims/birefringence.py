r"""Tier 2 (frontier): cosmic birefringence -- the parity-odd CMB fingerprint of
Eternal Dawn.

A parity-conserving sky (LambdaCDM) has exactly four CMB spectra: TT, EE, BB, TE.
The parity-ODD cross-spectra EB and TB vanish identically. Cosmic birefringence -- a
rotation of the linear-polarization plane by an angle beta on the way from last
scattering -- mixes E into B and lights them up:

    C_l^{EB} = (1/2) sin(4 beta) (C_l^{EE} - C_l^{BB}),
    C_l^{TB} =        sin(2 beta)  C_l^{TE}.

Because EE and TE are large, even a fraction-of-a-degree beta leaves a measurable
EB/TB -- which is how Planck and ACT have reported beta ~ 0.2-0.3 deg at 2.9-3.6
sigma. LambdaCDM predicts beta = 0 exactly, so any nonzero beta is parity violation
in the photon sector.

Why Eternal Dawn expects beta != 0. The framework carries a genuine global parity
violation (the matter/arrow choice co-sourced at the bounce; a chirally pure
lineage), and Einstein--Cartan gravity supplies a concrete agent: the totally
antisymmetric part of torsion is a pseudoscalar axial vector S_mu, which couples to
the photon like a Chern--Simons term and rotates polarization -- where LambdaCDM must
ADD an axion by hand, a relic axial-torsion background is already present. Two
components follow:

  * an ISOTROPIC beta_0 from a homogeneous parity-violating background -- a single
    sky-wide rotation, phenomenologically degenerate with the axion interpretation;

  * an ANISOTROPIC, AXIS-ALIGNED beta(n) from the inherited parent spin axis -- a
    low-multipole modulation locked to the SAME direction as the CMB "axis of evil"
    and the galaxy-spin handedness (Ch.6). This is the distinctive, uniquely-ED
    signature: not just beta != 0 (an axion gives that too) but a beta-anisotropy
    whose axis coincides with the other parity anomalies.

HONESTY ON THE MAGNITUDE. The framework does not yet predict the VALUE of beta: it is
set by the relic axial-torsion amplitude and its photon coupling, neither computed
(that is the owed calculation, Ch.8). What is computed here, from the real CAMB
EE/TE/BB transfer, is the OBSERVABLE the framework predicts given beta -- the EB and
TB spectra -- and the consistency of the measured beta ~ 0.3 deg with a small relic
parity violation. The decisive, ED-specific test is the axis alignment of the
anisotropic part, not the number itself.
"""

from __future__ import annotations

import numpy as np

# Published isotropic-birefringence measurements (deg, 1-sigma). LambdaCDM = 0.
#   (label, beta, sigma, citation-key)
BETA_MEASUREMENTS = [
    ("Planck PR3 (Minami--Komatsu 2020)", 0.35, 0.14, "minami2020"),
    ("Planck PR4 (Diego-Palazuelos 2022)", 0.30, 0.11, "diegopalazuelos2022"),
    ("Planck PR4 + WMAP (Eskilt--Komatsu 2022)", 0.34, 0.094, "eskilt2022"),
]


def _deg2rad(beta_deg):
    return np.deg2rad(np.asarray(beta_deg, dtype=float))


def eb_spectrum(ee, bb, beta_deg):
    """Birefringence-induced EB spectrum: C_l^{EB} = 1/2 sin(4 beta)(EE - BB).
    Inputs are D_l [microK^2]; returns D_l^{EB} [microK^2]."""
    b = _deg2rad(beta_deg)
    return 0.5 * np.sin(4.0 * b) * (np.asarray(ee) - np.asarray(bb))


def tb_spectrum(te, beta_deg):
    """Birefringence-induced TB spectrum: C_l^{TB} = sin(2 beta) TE."""
    b = _deg2rad(beta_deg)
    return np.sin(2.0 * b) * np.asarray(te)


def eb_peak_amplitude(ee, bb, beta_deg):
    """Peak |D_l^{EB}| [microK^2] for a given beta -- the headline EB amplitude."""
    return float(np.nanmax(np.abs(eb_spectrum(ee, bb, beta_deg))))


def combined_beta(measurements=None):
    """Inverse-variance combination of the isotropic-beta measurements.
    Returns (beta, sigma, significance) in degrees and sigma. NOTE: these analyses
    share Planck data and are correlated, so the true joint significance is lower
    than this naive combination; we report it only as an indicative band."""
    meas = measurements or BETA_MEASUREMENTS
    betas = np.array([m[1] for m in meas])
    sigs = np.array([m[2] for m in meas])
    w = 1.0 / sigs**2
    beta = float(np.sum(w * betas) / np.sum(w))
    sigma = float(np.sqrt(1.0 / np.sum(w)))
    return beta, sigma, beta / sigma


def axis_anisotropy_map(beta0_deg, beta1_deg, axis_lonlat=(260.0, 60.0), nside_deg=5.0):
    """A schematic axis-aligned birefringence field beta(n) = beta0 + beta1 * cos(theta),
    a dipole locked to the preferred axis (default: the CMB 'axis of evil' direction
    in Galactic (l,b)). Returns (lon, lat, beta_grid) on a regular grid [deg].

    This illustrates the ED-distinctive ANISOTROPIC component -- its dipole axis is
    predicted to coincide with the axis-of-evil / galaxy-spin axis (Ch.6). The
    amplitude beta1 is not derived; it is shown as a shape, not a prediction."""
    lon = np.arange(-180.0, 180.0 + nside_deg, nside_deg)
    lat = np.arange(-90.0, 90.0 + nside_deg, nside_deg)
    L, B = np.meshgrid(np.deg2rad(lon), np.deg2rad(lat))
    a_lon, a_lat = np.deg2rad(axis_lonlat[0]), np.deg2rad(axis_lonlat[1])
    # angular distance from the axis on the sphere
    cos_theta = (np.sin(B) * np.sin(a_lat)
                 + np.cos(B) * np.cos(a_lat) * np.cos(L - a_lon))
    beta = beta0_deg + beta1_deg * cos_theta
    return lon, lat, beta
