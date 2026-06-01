r"""The void as a maximum-entropy eigenstate, and OGUs as its suppressed fluctuations.

The thermodynamic question "what is the stable state of the void?" (Ch7: time is not
needed, only equilibrium) has a clean answer, and it ties three results together.

The eigenstate is EMPTY de Sitter space. Maximise the total horizon entropy of an
infinite self-gravitating vacuum with positive Lambda: in Schwarzschild-de Sitter the
generalised entropy is S = (A_bh + A_cosmo)/4 (Planck units), and inserting ANY mass
shrinks the cosmological horizon faster than it grows a black-hole horizon, so the
sum is MAXIMISED by the empty state. Empty de Sitter, with the Gibbons-Hawking entropy

    S_dS = pi / (l_p^2 H_Lambda^2)   ( = A_cosmo/4 ),

is therefore the maximum-entropy configuration -- the stable state, the eigenstate.

OGUs are NOT the eigenstate; they are its fluctuation spectrum. The probability of a
fluctuation that inserts an OGU is the Boltzmann/Einstein weight of its entropy
deficit, P ~ e^{Delta S} with Delta S < 0:

    Delta S(M) = S_SdS(M) - S_dS = -(A_cosmo deficit) ~ -2 pi M c^2/(hbar H_Lambda),

and that deficit is EXACTLY the OG-nucleation instanton action I (Ch5/Ch ?, the de
Sitter Boltzmann factor): Delta S = -I, so the nucleation rate lambda = e^{-I} is the
horizon-entropy suppression of the fluctuation. The thermodynamic eigenequation and
the nucleation instanton are one statement.

This unifies three results: (i) the supraverse is DILUTE because empty is
max-entropy, so every OGU is a rare entropy debt (I >> I_crit); (ii) the birth rate
beta carries the Boltzmann weight e^{-I} of that debt; (iii) each heat-dead interior
relaxes back to the same empty-de-Sitter eigenstate (low Weyl, max entropy) and
re-fluctuates -- the eigenstate is also the recursive attractor. It does NOT derive the
seed mass (set by rho_C / the bounce condition); it provides the variational principle
the equilibrium picture was implicitly using.
"""

from __future__ import annotations

import numpy as np

from . import constants as k
from . import sds


def hubble_lambda(Omega_Lambda: float = None, H0: float = None) -> float:
    OL = k.Omega_Lambda if Omega_Lambda is None else Omega_Lambda
    H = k.H0 if H0 is None else H0
    return H * np.sqrt(OL)


def de_sitter_horizon_entropy(H_lambda: float = None) -> float:
    """Gibbons-Hawking entropy of empty de Sitter, S/k_B = pi c^5/(hbar G H^2)
    = A_cosmo/(4 l_p^2). The maximum-entropy (eigenstate) value."""
    H = hubble_lambda() if H_lambda is None else H_lambda
    return np.pi * k.c**5 / (k.hbar * k.G * H**2)


def entropy_deficit(mu_frac: float) -> float:
    """Dimensionless horizon-entropy deficit S_dS - S_SdS, in units 1/(l_p^2 H^2),
    for a mass at fraction mu_frac of the Nariai mass. Positive: inserting mass
    LOWERS the total entropy, so empty de Sitter is the maximum."""
    mu = mu_frac * sds.MU_NARIAI
    h = sds.sds_horizons(mu)
    if h is None:
        return float("nan")
    rb, rc = h
    return float(np.pi * (1.0 - rb**2 - rc**2))   # S_empty(=pi) - S_SdS


def empty_is_maximum(n: int = 50) -> bool:
    """Verify empty de Sitter maximises the total horizon entropy: the deficit is
    positive (entropy drops) for every sub-Nariai mass."""
    fr = np.linspace(1e-4, 0.999, n)
    return bool(np.all([entropy_deficit(f) > 0 for f in fr]))


def deficit_small_mass(mu_frac: float) -> float:
    """Small-mass limit of the deficit, 2 pi mu (units 1/l_p^2H^2): equals the
    de Sitter Boltzmann exponent, i.e. the instanton action I in these units."""
    return 2.0 * np.pi * mu_frac * sds.MU_NARIAI


def nucleation_weight_is_entropy_deficit() -> bool:
    """The nucleation suppression e^{-I} equals e^{Delta S}: the instanton action is
    the horizon-entropy deficit. (Confirmed by construction; both = 2 pi M c^2/hbar H.)
    Here we check the small-mass deficit matches the exact deficit to leading order."""
    f = 1e-3
    return bool(abs(entropy_deficit(f) / deficit_small_mass(f) - 1.0) < 1e-2)
