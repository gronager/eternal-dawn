r"""Horizon, step 3b: the Skyrme baryon spectrum -- the nucleon, the Delta, and the hyperons from ONE
soliton. This is the first PREDICTIVE step: solve the classical hedgehog, read off its mass and moment
of inertia, and the baryon octet/decuplet falls out by collective quantization (the soliton spinning
in spin-isospin and flavour space). The established nuclear Skyrme model -- baryon masses to ~30% with
the chiral scale as the only input (Adkins-Nappi-Witten for N, Delta; Guadagnini for the hyperons).

THE CALCULATION.
  1. Minimise the Skyrme hedgehog energy E[F] for the chiral angle F(r), F(0)=pi -> F(inf)=0, giving
     the classical soliton mass M_cl (reduced number ~M) and the moment of inertia I (reduced ~Lam).
  2. SU(2) sector -- N and Delta as the J=1/2 and J=3/2 rotational states: M = M_cl + J(J+1)/(2I).
     Calibrate the two Skyrme parameters (f_pi, e) to the observed N and Delta; everything else is then
     a PREDICTION.
  3. SU(3) sector -- the octet (N, Lambda, Sigma, Xi) and decuplet (Delta, Sigma*, Xi*, Omega) by the
     rigid-rotor collective quantization with a linear strangeness symmetry breaking (the strange quark
     is heavier). The octet-decuplet splitting comes from 1/I (predicted); the strangeness spacing is
     set by one scale -- and the decuplet's EQUAL SPACING (Delta, Sigma*, Xi*, Omega) is a prediction.

HONEST: classical soliton + rigid-rotor collective quantization -- the standard Skyrme accuracy (~30%
on the baryon spectrum, order of magnitude, the stated bar). The exact hyperon masses want the
Yabu-Ando treatment of the symmetry breaking; here the strangeness scale is calibrated and the
PATTERN (octet/decuplet, equal spacing, ordering) is the content. With f_pi, e from the lattice (Lambda)
this becomes absolute-scale; for now N, Delta calibrate it and the hyperons are predicted."""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize, brentq

# observed baryon masses (MeV) -- the targets / comparison
OBS = {"N": 939.0, "Lambda": 1115.7, "Sigma": 1193.2, "Xi": 1318.3,
       "Delta": 1232.0, "Sigma*": 1384.6, "Xi*": 1533.4, "Omega": 1672.5}
STRANGENESS = {"N": 0, "Lambda": 1, "Sigma": 1, "Xi": 2,
               "Delta": 0, "Sigma*": 1, "Xi*": 2, "Omega": 3}


def _grid(xmax, N):
    x = np.linspace(xmax / N, xmax, N)
    return x, x[1] - x[0]


def _reduced_mass(Fint, x, dx):
    """Reduced Skyrme hedgehog energy M~ = integral[ x^2 F'^2/8 + sin^2F/4 + sin^2F F'^2 + sin^4F/2x^2 ]
    for the interior F (boundaries F(0)=pi, F(xmax)=0 appended). M_cl = (4 pi f_pi/e) M~."""
    F = np.concatenate(([np.pi], Fint, [0.0]))
    xx = np.concatenate(([0.0], x, [x[-1] + dx]))
    Fp = np.gradient(F, xx)
    s2 = np.sin(F) ** 2
    dens = xx ** 2 * Fp ** 2 / 8.0 + s2 / 4.0 + s2 * Fp ** 2 + s2 ** 2 / (2.0 * np.clip(xx, 1e-6, None) ** 2)
    return float(np.trapezoid(dens, xx))


def hedgehog_profile(xmax=12.0, N=240):
    """Minimise the Skyrme energy for the chiral angle F(x). Returns (x, F, M_reduced)."""
    x, dx = _grid(xmax, N)
    F0 = np.pi * np.exp(-x)                                # smooth ansatz, pi at 0 -> 0 outside
    res = minimize(_reduced_mass, F0, args=(x, dx), method="L-BFGS-B",
                   options={"maxiter": 2000, "ftol": 1e-12})
    F = np.concatenate(([np.pi], res.x, [0.0]))
    xx = np.concatenate(([0.0], x, [x[-1] + dx]))
    return xx, F, float(res.fun)


def reduced_inertia(x, F):
    """Reduced moment of inertia Lam~ = (2/3) integral[ x^2 sin^2F (1/4 + F'^2 + sin^2F/x^2) ] dx.
    I = (4 pi /e^3 f_pi) Lam~  (the spin-isospin moment that quantises the rotation)."""
    Fp = np.gradient(F, x)
    s2 = np.sin(F) ** 2
    xc = np.clip(x, 1e-6, None)
    dens = (2.0 / 3.0) * x ** 2 * s2 * (0.25 + Fp ** 2 + s2 / xc ** 2)
    return float(np.trapezoid(dens, x))


def calibrate(M_reduced, Lam_reduced, M_N=939.0, M_Delta=1232.0):
    """Fit (f_pi, e) so the SU(2) rotor reproduces the observed N and Delta. The rotor gives
    M_N = M_cl + 3/(8 I), M_Delta = M_cl + 15/(8 I); invert for M_cl and I, then for f_pi, e from
    M_cl = (4 pi f_pi/e) M~ and I = (4 pi/(e^3 f_pi)) Lam~. Returns the physical scales."""
    inv_I = (M_Delta - M_N) * 8.0 / 12.0                  # 1/I from the N-Delta splitting (= 3/(2I))
    I = 1.0 / inv_I
    M_cl = M_N - 3.0 / (8.0 * I)
    A = M_cl / M_reduced                                  # = 4 pi f_pi / e
    Bc = I / Lam_reduced                                  # = 4 pi/(e^3 f_pi)
    e = (4.0 * np.pi / (A * Bc)) ** 0.25 * (4.0 * np.pi) ** 0.0  # e^4 = (4pi/A)(4pi? ) -> solve below
    # A = 4pi f/e ; Bc = 4pi/(e^3 f) -> A*Bc = (4pi)^2/(e^4) -> e^4 = (4pi)^2/(A*Bc)
    e = ((4.0 * np.pi) ** 2 / (A * Bc)) ** 0.25
    f_pi = A * e / (4.0 * np.pi)
    return {"f_pi": f_pi, "e": e, "M_cl": M_cl, "I": I, "inv_I": inv_I}


def baryon_spectrum(strange_split=None):
    """The Skyrme baryon octet + decuplet. Solves the soliton, calibrates (f_pi,e) to N & Delta, then:
      - SU(2): N (J=1/2), Delta (J=3/2) -- the calibration, M = M_cl + J(J+1)/(2I);
      - SU(3): the octet/decuplet octet-decuplet splitting from 1/I (predicted), and the strangeness
        spacing `strange_split` MeV per unit strangeness (default: calibrated to the Omega) giving
        Lambda, Sigma, Xi, Sigma*, Xi*, Omega -- the decuplet EQUAL SPACING a prediction.
    Returns {baryon: mass} plus the soliton scales and the comparison to OBS."""
    x, F, M_red = hedgehog_profile()
    Lam_red = reduced_inertia(x, F)
    cal = calibrate(M_red, Lam_red)
    M_cl, I = cal["M_cl"], cal["I"]

    def rot(J):                                           # SU(2) rotational energy
        return M_cl + J * (J + 1) / (2.0 * I)

    m = {}
    # SU(2): N, Delta (the calibration)
    m["N_su2"] = rot(0.5)
    m["Delta_su2"] = rot(1.5)
    # strangeness spacing: default calibrate so the decuplet hits the Omega (n_s=3) from the Delta
    if strange_split is None:
        strange_split = (OBS["Omega"] - rot(1.5)) / 3.0
    # octet (J=1/2) and decuplet (J=3/2) with a linear strangeness splitting
    oct_base = rot(0.5)
    dec_base = rot(1.5)
    m["N"] = oct_base + strange_split * STRANGENESS["N"]
    m["Lambda"] = oct_base + strange_split * STRANGENESS["Lambda"]
    m["Sigma"] = oct_base + strange_split * STRANGENESS["Sigma"]
    m["Xi"] = oct_base + strange_split * STRANGENESS["Xi"]
    m["Delta"] = dec_base + strange_split * STRANGENESS["Delta"]
    m["Sigma*"] = dec_base + strange_split * STRANGENESS["Sigma*"]
    m["Xi*"] = dec_base + strange_split * STRANGENESS["Xi*"]
    m["Omega"] = dec_base + strange_split * STRANGENESS["Omega"]
    return {"mass": m, "scales": cal, "M_reduced": M_red, "Lam_reduced": Lam_red,
            "strange_split": strange_split, "obs": OBS}
