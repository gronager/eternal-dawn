r"""1-D string tension: does the gravity-torsion vacuum confine? (qualitative)

The strong sector's deepest open question (Part II's downgrade): the overlap/exchange
forces all SCREEN (-> 0), but real quarks are CONFINED (V ~ sigma r, never zero).
Confinement is a non-perturbative FLUX TUBE, not a pairwise overlap. This module tests
the standard mechanism by which it could arise here -- the DUAL SUPERCONDUCTOR: if the
chiral condensate (the sigma field the soliton forms) expels colour-electric flux, the
flux is squeezed into a Nielsen-Olesen vortex (a 1-D tube), whose energy is proportional
to its length -> confinement, with a computable string tension.

The vortex (Abelian-Higgs, winding n=1, units v=1, m_W^2=2):
  f(r) -- condensate amplitude: 0 on the axis (normal core), 1 outside (superconductor);
  a(r) -- gauge profile: the colour-electric flux, confined to the tube.
Energy/length sigma = 2pi integral[ r f'^2 + (1-a)^2 f^2/r + a'^2/(2r) + (beta/4) r (f^2-1)^2 ].
At the BPS point (beta=2, m_H=m_W) the topological bound is saturated: sigma = 2 pi v^2
EXACTLY (the code's validation). For general beta, sigma ~ O(1) x 2pi v^2.

THE RESULT. A confining flux tube EXISTS with tension sigma proportional to the
condensate scale v^2 = f_pi^2 -- so the SAME condensate that makes mass (configurational,
Part II) also confines. Confinement gives V(L) = sigma L (rises forever), the linear
potential the overlap picture could not produce.

HONEST STATUS. This demonstrates the MECHANISM and the tension SCALE, conditional on the
key assumption: that the gravity-torsion vacuum is a dual superconductor (its condensate
expels colour flux via a dual Meissner effect). Whether it actually is -- monopole
condensation in the non-abelian connection -- is the genuinely non-perturbative question,
and it is exactly what a real LATTICE gauge-theory computation would decide (measure the
Wilson-loop area law / the dual order parameter for the Part I connection). So: confinement
is POSSIBLE and its scale is set by the condensate; whether the GT vacuum realises it is
owed, and is a clean, well-posed target for lattice. Qualitative, with a lattice path.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_bvp

_trapz = np.trapezoid


def solve_vortex(beta=2.0, R=14.0, N=3000):
    """Solve the Nielsen-Olesen vortex ODEs. Returns dict(r, f, a, energy_density,
    sigma, success). sigma is the string tension (energy per unit length)."""
    r = np.linspace(1e-4, R, N)

    def ode(r, y):
        f, fp, a, ap = y
        fpp = -fp / r + (1 - a) ** 2 * f / r**2 + 0.5 * beta * (f**2 - 1) * f
        app = ap / r - 2 * (1 - a) * f**2
        return np.vstack([fp, fpp, ap, app])

    def bc(ya, yb):
        return np.array([ya[0], ya[2], yb[0] - 1, yb[2] - 1])

    y0 = np.zeros((4, N))
    y0[0] = np.tanh(r)
    y0[2] = np.tanh(r / 2) ** 2
    sol = solve_bvp(ode, bc, r, y0, max_nodes=300000, tol=1e-7)
    f, fp, a, ap = sol.sol(r)
    dens = r * fp**2 + (1 - a) ** 2 * f**2 / r + ap**2 / (2 * r) + (beta / 4) * r * (f**2 - 1) ** 2
    sigma = 2 * np.pi * float(_trapz(dens, r))
    # the magnetic (colour-electric) field B(r) ~ a'/r, localised in the tube
    B = ap / r
    return {"r": r, "f": f, "a": a, "B": B, "energy_density": 2 * np.pi * dens,
            "sigma": sigma, "success": bool(sol.success), "beta": beta}


def string_tension(beta=2.0):
    """String tension sigma (units of v^2) for Ginzburg-Landau parameter beta."""
    return solve_vortex(beta)["sigma"]


def bps_tension():
    """At the BPS point (beta=2 here, m_H=m_W) the tension saturates: sigma = 2 pi v^2."""
    return string_tension(2.0)


def confining_potential(L, sigma):
    """V(L) = sigma * L: the linear confining potential between two sources joined by a
    flux tube (rises forever -> quarks cannot be isolated)."""
    return sigma * np.asarray(L, dtype=float)
