r"""Pinning the OG birth rate beta -- as the de Sitter nucleation rate.

The recursive cycle (cycle.py) reinterprets the OG birth rate: OGUs nucleate in the
heat-dead de Sitter void of an older universe. So beta is a de Sitter vacuum
nucleation rate, and its scale is set -- not free -- by the void's dark energy.
Writing it as a dimensionless rate per de Sitter horizon 4-volume,

    beta = lambda * H_Lambda^4 / c^3,

the PREFACTOR H_Lambda^4/c^3 is fixed by the observed dark energy (by typicality the
void's Lambda is ~ ours), H_Lambda = H_0 sqrt(Omega_Lambda) ~ 1.8e-18 s^-1, giving

    H_Lambda^4 / c^3 ~ 4e-97 m^-3 s^-1   (one "attempt" per Hubble 4-volume).

The one remaining unknown is the dimensionless lambda = exp(-I), with I the action
of the rho_C-crossing gravitational instanton (quantum gravity at the bounce, Q4).
But lambda has a SHARP critical value: de Sitter bubble percolation (Guth-Weinberg)
sets lambda_crit ~ 0.24, equivalently an instanton action I_crit = -ln(lambda_crit)
~ 1.4. This decides the entire CHARACTER of the supraverse:

  * lambda > lambda_crit (I < 1.4): OGUs nucleate faster than de Sitter dilutes ->
    they percolate -> the packed polygonal foam that coarsens (void_foam.py).
  * lambda < lambda_crit (I > 1.4): the void inflates between nucleations -> OGUs
    stay isolated, round, never impinging -> a dilute scatter of island universes
    (the eternal-inflation regime); no coarsening, no merging.

So beta is pinned in scale (the dark-energy prefactor) and its decisive structure
(lambda_crit ~ 0.24, I_crit ~ 1.4) is identified; the regime hinges on the bounce
instanton I, still open. Our being a BHU1 inside an OGU is compatible with either
regime (the internal black-hole channel is independent of how OGUs are spaced).
"""

from __future__ import annotations

import numpy as np

from . import constants as k

LAMBDA_CRIT = 0.24                 # de Sitter percolation threshold (Guth-Weinberg)


def hubble_lambda(Omega_Lambda: float = None, H0: float = None) -> float:
    """de Sitter expansion rate from dark energy, H_Lambda = H_0 sqrt(Omega_Lambda)
    [s^-1]; by typicality the void's value is ~ ours."""
    OL = k.Omega_Lambda if Omega_Lambda is None else Omega_Lambda
    H = k.H0 if H0 is None else H0
    return H * np.sqrt(OL)


def nucleation_prefactor(**kw) -> float:
    """H_Lambda^4 / c^3 [m^-3 s^-1]: the de Sitter attempt rate per Hubble
    4-volume -- the scale of beta, fixed by the observed dark energy."""
    return hubble_lambda(**kw) ** 4 / k.c**3


def beta(lam: float, **kw) -> float:
    """OG birth rate beta = lambda * H_Lambda^4/c^3 [m^-3 s^-1]."""
    return lam * nucleation_prefactor(**kw)


def instanton_action(lam: float) -> float:
    """The bounce instanton action implied by a dimensionless rate, I = -ln(lambda)."""
    return -np.log(lam)


def percolates(lam: float) -> bool:
    """True if OGUs nucleate fast enough to percolate into a packed foam."""
    return lam >= LAMBDA_CRIT


def regime(lam: float) -> str:
    return ("percolating: packed polygonal foam (coarsens)" if percolates(lam)
            else "dilute: isolated round island universes (eternal-inflation-like)")


def ogu_separation_hubble(lam: float) -> float:
    """Mean comoving OGU separation in de Sitter horizon radii, ~ lambda^{-1/4}
    (the 4-volume nucleation spacing). ~1 at percolation, vast for tiny lambda."""
    return lam ** (-0.25)
