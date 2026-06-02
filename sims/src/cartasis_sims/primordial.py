r"""Tier 2 (first result): the primordial spectrum from the Einstein--Cartan bounce.

The decisive test of the framework (Ch.9): does the bounce produce a near-scale-
invariant primordial spectrum, n_s ~ 0.96, without inflation? This module evolves
perturbation modes through the bounce and reads off the spectral tilt.

Mechanism (the "matter bounce"). Perturbation modes start sub-horizon deep in the
contracting phase, exit the comoving Hubble horizon during contraction (the bounce
analogue of inflation's horizon exit), freeze, cross the bounce, and re-enter in
expansion. The frozen super-horizon amplitude sets the primordial spectrum. For a
contracting phase with constant equation of state w, the conformal scale factor is
a power law a ~ |eta|^p with

    p = 2/(1 + 3w),

and the tensor (gravitational-wave) spectrum is P_T(k) ~ k^{n_T} with

    n_T = 3 - |2p - 1|.

Matter-dominated contraction (w = 0, p = 2) gives n_T = 0 -- EXACTLY scale-
invariant -- which is the observed near-scale-invariance, reproduced with no
inflaton. Radiation (w = 1/3, p = 1) gives n_T = 2 (steeply blue). So the bounce
explains the scale-invariance IF the contraction is matter-dominated; the observed
slight red tilt (n_s = 0.965) corresponds to a small positive effective w.

What the Einstein--Cartan bounce adds is the NON-SINGULAR bounce itself: the matter
bounce scenario otherwise needs exotic ghost matter to turn the contraction
around, whereas EC torsion does it automatically at rho_C (Ch.2). The scalar
spectrum follows the same mechanism; here we compute the cleanly-defined tensor
spectrum and validate it against the analytic tilt. The precise n_s (the small red
tilt) and the tensor-to-scalar ratio r need the detailed near-matter contraction
history -- the remaining Tier-2 work.

Modes are evolved with the Mukhanov--Sasaki equation
    mu_k'' + (k^2 - a''/a) mu_k = 0,   h_k = mu_k / a,
in conformal time, with Bunch--Davies initial data deep in contraction, and the
spectrum read while super-horizon (at the bounce), P_T(k) ~ k^3 |h_k|^2.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

ETA_B = 1.0           # bounce conformal duration (units); sets a_min = 1
N_S_OBS = 0.9649      # Planck 2018 scalar spectral index


def p_from_w(w: float) -> float:
    """Conformal power-law index of the contracting scale factor a ~ |eta|^p."""
    return 2.0 / (1.0 + 3.0 * w)


def analytic_tilt(p: float) -> float:
    """Tensor spectral index of the matter-bounce family: n_T = 3 - |2p - 1|."""
    return 3.0 - abs(2.0 * p - 1.0)


def scale_factor(eta, p: float):
    """Bounce background a(eta) = (1 + (eta/eta_b)^2)^{p/2}; a -> |eta|^p far out,
    smooth non-singular bounce a_min = 1 at eta = 0 (the Einstein--Cartan turn)."""
    return (1.0 + (np.asarray(eta) / ETA_B) ** 2) ** (p / 2.0)


def a_second_over_a(eta, p: float):
    """a''/a in conformal time for the bounce family (analytic)."""
    u = 1.0 + (eta / ETA_B) ** 2
    a = u ** (p / 2.0)
    app = (p / ETA_B**2 * u ** (p / 2.0 - 1.0)
           + p * eta / ETA_B**2 * (p / 2.0 - 1.0) * u ** (p / 2.0 - 2.0)
           * 2.0 * eta / ETA_B**2)
    return app / a


def mode_power(k: float, p: float, eta0: float = -400.0) -> float:
    """|h_k|^2 read super-horizon (at the bounce, eta=0) for one mode k.

    Bunch--Davies sub-horizon initial data deep in contraction; evolve the
    Mukhanov--Sasaki equation to the bounce where the mode is frozen."""
    def rhs(eta, y):
        mr, mi, dr, di = y
        A = a_second_over_a(eta, p)
        return [dr, di, -(k * k - A) * mr, -(k * k - A) * mi]

    phase = np.exp(-1j * k * eta0) / np.sqrt(2.0 * k)
    y0 = [phase.real, phase.imag, (-1j * k * phase).real, (-1j * k * phase).imag]
    sol = solve_ivp(rhs, (eta0, 0.0), y0, rtol=1e-9, atol=1e-12, max_step=0.3)
    mr, mi = sol.y[0][-1], sol.y[1][-1]
    return (mr * mr + mi * mi) / scale_factor(0.0, p) ** 2


def power_spectrum(ks, p: float):
    """Dimensionless tensor power P_T(k) ~ k^3 |h_k|^2 over an array of modes."""
    ks = np.asarray(ks, dtype=float)
    return np.array([k**3 * mode_power(k, p) for k in ks])


def tilt_numeric(p: float, ks=None) -> float:
    """Numerical tensor tilt n_T = d ln P_T/d ln k, fit over super-horizon modes."""
    if ks is None:
        ks = np.logspace(-2.0, -0.7, 9)      # all super-horizon at the bounce
    P = power_spectrum(ks, p)
    return float(np.polyfit(np.log(ks), np.log(P), 1)[0])
