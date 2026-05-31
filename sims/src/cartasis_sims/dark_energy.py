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
# Phantom-crossing dark energy from a peaked parent-injection history
# ---------------------------------------------------------------------------
# Dark energy here is sourced by matter injected through the persistent Cartasis
# membrane from the accreting parent (Ch. 6). If we absorb the source into an
# *effective* equation of state (no explicit source term), the continuity
# equation rho_dot + 3H(1+w_eff) rho = 0 gives the exact identity
#
#     w_eff(a) = -1 - (1/3) d ln rho_DE / d ln a.
#
# So rho_DE *rising* with a (injection winning over dilution) => w_eff < -1
# (phantom); rho_DE *falling* (dilution winning) => w_eff > -1; rho_DE flat =>
# w_eff = -1. A parent whose net injection peaked at scale factor a_p therefore
# makes rho_DE(a) a hump that crosses from phantom (past, a<a_p) to
# quintessence-like (now, a>a_p), crossing w=-1 exactly at a_p. This is the DESI
# pattern, derived rather than fitted: the crossing is a *consequence* of
# injection that rose and then saturated, not a free sign choice.

def rho_de_injection(a, a_p=0.75, beta=1.3, rho0=1.0):
    r"""Peaked (log-normal in a) dark-energy density from parent injection.

    rho_DE(a) = rho0 * exp[-beta (ln a - ln a_p)^2], peaked at a=a_p, where a_p
    is the scale factor at which the parent's net injection (accretion minus the
    baby's expansion dilution) turned over.
    """
    a = np.asarray(a, dtype=float)
    return rho0 * np.exp(-beta * (np.log(a) - np.log(a_p)) ** 2)


def w_eff_injection(a, a_p=0.75, beta=1.3):
    r"""Effective w(a) for the peaked-injection density, from the exact identity
    w_eff = -1 - (1/3) d ln rho/d ln a. For the log-normal bump
    d ln rho/d ln a = -2 beta (ln a - ln a_p), so

        w_eff(a) = -1 + (2 beta / 3) (ln a - ln a_p),

    which is < -1 for a < a_p (phantom past) and > -1 for a > a_p (today).
    """
    a = np.asarray(a, dtype=float)
    return -1.0 + (2.0 * beta / 3.0) * (np.log(a) - np.log(a_p))


def injection_from_cpl(w0, wa):
    r"""Invert (w0, wa) -> (a_p, beta) for the peaked-injection model.

    Linearizing w_eff about a=1 (ln a ~ -(1-a)) gives CPL with
        w0 = -1 - (2 beta/3) ln a_p,   wa = -2 beta/3.
    Hence beta = -3 wa / 2 and ln a_p = (w0 + 1)/wa, so a_p = exp((w0+1)/wa).
    Returns a_p, beta, and the phantom-crossing redshift z_cross = 1/a_p - 1.
    """
    beta = -1.5 * wa
    a_p = float(np.exp((w0 + 1.0) / wa)) if wa != 0 else float("nan")
    z_cross = 1.0 / a_p - 1.0 if a_p == a_p else float("nan")
    return {"a_p": a_p, "beta": beta, "z_cross": z_cross}


def logistic_accretion(a, a_mid=0.4, k=6.0, m_inf=1.0):
    r"""Parent mass fraction vs the baby's scale factor: Eddington-like (near
    exponential) growth that saturates as fuel depletes -- a logistic in ln a,
    M_p(a)/M_inf = 1/(1 + (a/a_mid)^{-k}). Its log-derivative (fractional
    accretion rate) peaks near a_mid and decays; that turnover is what seeds the
    injection bump and hence the phantom->quintessence crossing.
    """
    a = np.asarray(a, dtype=float)
    return m_inf / (1.0 + (a / a_mid) ** (-k))


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
