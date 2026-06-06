r"""L2 -- the string tension from the condensate: sigma = 2 pi v^2.

The book's confinement section (sec:confinement) predicts that the SAME chiral condensate v
that sets configurational mass also sets the confinement scale: a Nielsen--Olesen (Abrikosov)
vortex of that condensate squeezes colour flux into a string of tension

    sigma = 2 pi v^2          (winding-1, BPS / Bogomolny bound).

"One substrate does both." This module is the L2 core:

  * bps_string_tension(v)      -- the prediction sigma = 2 pi v^2 from the chiral-soliton v.
  * ano_vortex_tension(beta)   -- the ACTUAL vortex tension (in units of 2 pi v^2) as a function
                                  of the Ginzburg--Landau parameter beta = (m_scalar/m_vector)^2,
                                  by solving the vortex BVP. At BPS (beta=1) it saturates the
                                  bound exactly (returns 1.0); beta<1 (type I) gives <1, beta>1
                                  (type II) gives >1. So 2 pi v^2 is the prediction AT BPS, and
                                  whether the gravity-torsion vacuum sits there (with beta from
                                  the propagating-torsion coupling, lattice target L5) is what
                                  makes the relation exact or merely close.
  * sigma_fpi_ratio(sigma, fpi)-- the dimensionless, scale-free LATTICE test, sigma / f_pi^2,
                                  predicted = 2 pi. The GH200 measures sigma (static potential,
                                  measure_potential) AND f_pi (axial correlator) on the SAME
                                  chiral-breaking configs and forms this ratio. A pure-gauge
                                  (quenched) sigma cannot enter -- it has no f_pi.

NOT used here: the quenched pure-gauge sigma from out/puregauge (no f_pi, no common scale).
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_bvp

TWO_PI = 2.0 * np.pi


def bps_string_tension(v):
    """sigma = 2 pi v^2: the winding-1 BPS flux-tube tension, the framework's L2 prediction.
    v is the chiral condensate amplitude (f_pi) from chiral_soliton.solve_chiral()."""
    return TWO_PI * np.asarray(v, dtype=float) ** 2


def sigma_fpi_ratio(sigma, fpi):
    """The dimensionless L2 observable sigma / f_pi^2, predicted to equal 2 pi. Measure both on
    the SAME chiral-breaking lattice configs (sigma from the static potential, f_pi from the
    axial correlator) -- this is the scale-free, make-or-break number."""
    return float(np.asarray(sigma, dtype=float) / np.asarray(fpi, dtype=float) ** 2)


def ano_vortex_tension(beta, x_max=60.0, n=6000, max_nodes=200000):
    """Tension of a winding-1 Abrikosov--Nielsen--Olesen vortex of the condensate, in units of
    2 pi v^2, vs the Ginzburg--Landau parameter beta = (m_scalar/m_vector)^2. Dimensionless BVP
    (radius x = m_vector * r):

        f'' + f'/x - (1-a)^2 f / x^2 + (beta/2) f (1 - f^2) = 0      (scalar profile f: 0 -> 1)
        a'' - a'/x + (1-a) f^2 = 0                                    (gauge profile a: 0 -> 1)

    tension / (2 pi v^2) = \int_0^inf x dx [ f'^2 + (1-a)^2 f^2/x^2 + a'^2/x^2 + (beta/4)(1-f^2)^2 ].

    At BPS (beta=1) this saturates to exactly 1 (sigma = 2 pi v^2)."""
    x = np.linspace(1.0e-3, x_max, n)
    f0 = np.tanh(x)
    a0 = x ** 2 / (1.0 + x ** 2)
    y0 = np.vstack([f0, np.gradient(f0, x), a0, np.gradient(a0, x)])

    def odes(x, y):
        f, fp, a, ap = y
        fpp = -fp / x + (1.0 - a) ** 2 * f / x ** 2 - (beta / 2.0) * f * (1.0 - f ** 2)
        app = ap / x - (1.0 - a) * f ** 2
        return np.vstack([fp, fpp, ap, app])

    def bc(yl, yr):
        return np.array([yl[0], yl[2], yr[0] - 1.0, yr[2] - 1.0])

    sol = solve_bvp(odes, bc, x, y0, max_nodes=max_nodes, tol=1e-7)
    if not sol.success:
        raise RuntimeError(f"ANO vortex BVP did not converge (beta={beta}): {sol.message}")
    f, fp, a, ap = sol.sol(x)
    integrand = x * (fp ** 2 + (1.0 - a) ** 2 * f ** 2 / x ** 2
                     + ap ** 2 / x ** 2 + (beta / 4.0) * (1.0 - f ** 2) ** 2)
    trapezoid = getattr(np, "trapezoid", getattr(np, "trapz", None))
    return float(trapezoid(integrand, x))


def predicted_sigma_from_soliton(solve_chiral, **kwargs):
    """Convenience: run the chiral soliton, read off f_pi = v, and return the L2 prediction
    sigma = 2 pi v^2 together with v. `solve_chiral` is cartasis_sims.chiral_soliton.solve_chiral."""
    out = solve_chiral(**kwargs)
    v = float(out["f_pi"]) if isinstance(out, dict) and "f_pi" in out else float(out)
    return {"v": v, "sigma_bps": float(bps_string_tension(v))}
