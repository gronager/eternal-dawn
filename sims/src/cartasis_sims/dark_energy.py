r"""Dark sector vs data: parent-accretion w(z), dark/baryon ratio, BH growth.

The framework's claim (Chapter 6): inside a baby universe, dark energy is the
response to the *parent* black hole's growth, so the effective equation of state
is not -1 but tracks the parent's accretion history. We build the simplest such
toy model and compare its (w0, wa) to the DESI DR2 CPL fits.

Toy model. Let the parent mass grow as a power of the baby's scale factor,
M_parent(a) ~ a^p with p > 0 while it is still accreting (p -> 0 as it saturates).
The membrane-driven expansion behaves like a dark-energy component whose pressure
softens as the parent's *fractional* growth rate d ln M / d ln a = p declines as
the baby expands. A minimal, dimensionless ansatz reproducing this is

    w(a) = -1 + eps * f_grow(a),    f_grow(a) = (a/a)*p_eff(a),

with p_eff decreasing toward the future (growth saturates). Expanded around a=1
in CPL form w(a) = w0 + wa (1 - a):

    w0 = -1 + eps*p0,      wa = +eps * dp/dln a |_0  (negative if growth saturates).

So a still-growing-but-decelerating parent gives w0 > -1 and wa < 0 -- the DESI
sign pattern -- without a cosmological constant. eps absorbs the (unknown)
membrane coupling; p0 and its log-derivative are the parent-growth inputs.
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Toy parent-accretion equation of state
# ---------------------------------------------------------------------------
def w_of_a(a, eps=0.5, p0=0.6, q=1.2):
    r"""Effective dark-energy w(a) from a saturating parent-growth rate.

    Parent fractional growth rate p_eff(a) = p0 * a^q with q>0 means growth was
    faster (per e-fold) in the past and saturates toward the future. Then
    w(a) = -1 + eps * p_eff(a).  At a=1: w0 = -1 + eps p0; dw/da = eps p0 q,
    so wa = -dw/da|_1 ... (CPL: w = w0 + wa(1-a) => wa = -dw/da|_1).
    """
    a = np.asarray(a, dtype=float)
    return -1.0 + eps * p0 * a**q


def cpl_from_toy(eps=0.5, p0=0.6, q=1.2):
    """Return (w0, wa) of the CPL match to the toy model at a=1."""
    w0 = -1.0 + eps * p0
    dw_da = eps * p0 * q                 # d w / d a at a=1
    wa = -dw_da                          # CPL: w(a)=w0+wa(1-a) -> wa=-dw/da|_1
    return w0, wa


def w_cpl(a, w0, wa):
    """Standard CPL equation of state w(a) = w0 + wa (1 - a)."""
    a = np.asarray(a, dtype=float)
    return w0 + wa * (1.0 - a)


def fit_toy_to(w0_target, wa_target, p0=0.6, q=None):
    """Pick (eps, q) so the toy model reproduces a target (w0, wa).

    From w0 = -1 + eps p0  =>  eps = (w0+1)/p0.
    From wa = -eps p0 q    =>  q  = -wa/(eps p0) = -wa/(w0+1).
    """
    eps = (w0_target + 1.0) / p0
    q_fit = -wa_target / (w0_target + 1.0) if (w0_target + 1.0) != 0 else np.nan
    return {"eps": eps, "p0": p0, "q": q_fit}


# ---------------------------------------------------------------------------
# Dark/baryon ratio and the filter
# ---------------------------------------------------------------------------
def filter_fraction_from_ratio(dm_to_baryon):
    """f = 1/(1+ratio): pass-fraction if DM is filter-rejected counterpart mass."""
    return 1.0 / (1.0 + dm_to_baryon)


# ---------------------------------------------------------------------------
# Black-hole growth: critical mass for accretion-beats-Hawking
# ---------------------------------------------------------------------------
def m_crit_growth(T_bath, G=6.674e-11, c=2.998e8, hbar=1.055e-34, kB=1.381e-23):
    """M above which T_Hawking < T_bath, so the hole grows: hbar c^3/(8 pi G kB T)."""
    import math
    return hbar * c**3 / (8.0 * math.pi * G * kB * T_bath)
