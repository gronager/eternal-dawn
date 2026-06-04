r"""Degeneracy pressure in the soliton: the self-bound relativistic Fermi ball.

The first soliton pass (soliton.py) put ONE fermion in a fixed well -- no Pauli
exclusion, no degeneracy pressure. But a real soliton is many fermions filling the
levels, and their degeneracy pressure is a genuine OUTWARD force (the same one that
holds up white dwarfs and neutron stars), countering the inward binding alongside
the torsion wall. This module adds it at the Thomas--Fermi (equation-of-state) level
and shows which way it moves the soliton.

Energy per particle (units m=1, hbar=c=1, g internal states):

    e(n) = e_kin(n)  -  a n  +  b n^2 ,

  * e_kin(n) = 3 f(x)/x^3,  x = k_F/m = (6 pi^2 n/g)^{1/3},
        f(x) = [x(2x^2+1)sqrt(1+x^2) - asinh(x)]/8   (SymPy-verified)
    -- the relativistic Fermi-gas energy: rest mass + DEGENERACY kinetic energy.
  * - a n : the attractive binding (deepens with density);
  * + b n^2 : the torsion (Hehl--Datta) rho^2 wall -- repulsive at high density.

Minimising e(n) gives the self-bound saturation density n0 and the binding per
particle. Toggling the degeneracy term off (e_kin -> 1, rest mass only) isolates the
Pauli effect. The finding: degeneracy LOWERS n0 and SHALLOWS the binding -- the
soliton becomes larger, less dense, and more uniform (a saturated 'drop' with a flat
interior and a surface), i.e. more box-like / scale-flat in its core. That is the
direction that helps the electroweak S parameter (a flatter, more conformal interior
is the 'walking' regime), so Pauli pushes the right way -- it does not by itself
deliver S < 0.1, but it is in the favourable column.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar


def fermi_f(x):
    """f(x) = [x(2x^2+1)sqrt(1+x^2) - asinh(x)]/8  (SymPy-verified)."""
    x = np.asarray(x, dtype=float)
    return (x * (2 * x**2 + 1) * np.sqrt(1 + x**2) - np.arcsinh(x)) / 8.0


def kF_of_n(n, g=2.0):
    """Fermi momentum k_F (units m=1) for number density n, g internal states."""
    return (6.0 * np.pi**2 * n / g) ** (1.0 / 3.0)


def e_kin_per_particle(n, g=2.0):
    """Relativistic Fermi-gas energy per particle (rest + degeneracy KE), units m."""
    n = np.asarray(n, dtype=float)
    x = kF_of_n(np.maximum(n, 1e-300), g)
    out = np.where(n > 0, 3.0 * fermi_f(x) / np.maximum(x, 1e-300) ** 3, 1.0)
    return out


def energy_per_particle(n, a=5.0, b=4.0, g=2.0, degeneracy=True):
    """e(n) = e_kin(n) - a n + b n^2.  degeneracy=False uses rest mass only (no Pauli
    KE) to isolate the effect of the exclusion pressure."""
    ekin = e_kin_per_particle(n, g) if degeneracy else 1.0
    return ekin - a * n + b * n**2


def saturation(a=5.0, b=4.0, g=2.0, degeneracy=True, n_hi=0.6):
    """Self-bound saturation: the density n0 minimising e(n), and the binding per
    particle B = m - e(n0) (positive if bound). Returns dict(n0, e0, binding)."""
    res = minimize_scalar(lambda n: energy_per_particle(max(n, 1e-9), a, b, g, degeneracy),
                          bounds=(1e-4, n_hi), method="bounded")
    n0 = float(res.x)
    e0 = float(energy_per_particle(n0, a, b, g, degeneracy))
    return {"n0": n0, "e0": e0, "binding": 1.0 - e0, "bound": e0 < 1.0}


def degeneracy_shift(a=5.0, b=4.0, g=2.0):
    """How Pauli degeneracy moves the soliton: ratios of (density, binding) with vs
    without the exclusion KE. n0 drops and binding shallows -> larger, flatter drop."""
    on = saturation(a, b, g, degeneracy=True)
    off = saturation(a, b, g, degeneracy=False)
    return {
        "n0_with": on["n0"], "n0_without": off["n0"],
        "n0_ratio": on["n0"] / off["n0"],
        "binding_with": on["binding"], "binding_without": off["binding"],
        "binding_ratio": on["binding"] / max(off["binding"], 1e-9),
    }
