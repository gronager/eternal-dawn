r"""The electroweak S parameter from the soliton sector -- leading order, honestly.

S is Part III's make-or-break. A rigorous value for a strongly-coupled sector is a
lattice-scale question (the field has not settled it for technicolor in forty years),
so this module does the honest, computable part and frames the rest precisely.

WHAT IS COMPUTED. The leading contribution is the constituent-fermion loop: the
vector minus axial-vector current correlator of a fermion that got its mass M = g v
from the condensate (chiral_soliton.py). The free-fermion spectral functions are

    rho_V(s) = (N_c/12 pi)(1 + 2 M^2/s) sqrt(1 - 4 M^2/s),
    rho_A(s) = (N_c/12 pi)(1 - 4 M^2/s) sqrt(1 - 4 M^2/s),     s > 4 M^2,

so the chiral-odd DIFFERENCE is positive and UV-finite,

    rho_V(s) - rho_A(s) = (N_c/12 pi)(6 M^2/s) sqrt(1 - 4 M^2/s) > 0,

which is why S > 0: the vector channel outweighs the axial. The leading-order value
of a chirally-broken doublet is the standard result

    S_leading = N_c / (6 pi)  per electroweak doublet,

(Peskin & Takeuchi 1992) -- about 0.16 for N_c = 3, and the colour/doublet count
multiplies it. That is the technicolor 'graveyard': above the LEP bound ~0.1.

WHAT IS NOT. The reduction below 0.1 needs the SUBLEADING resonance/meson-loop
(walking) corrections -- the V and A composite resonances rearranging the spectral
functions so the difference cancels (the Weinberg-sum-rule form,
S = 4 pi (f_pi/M_V)^2 (1 + M_V^2/M_A^2)). Whether the anharmonic, minimal,
Pauli-flattened soliton sector walks enough is a genuine open question requiring RPA
correlators -- beyond this module. So: leading S is in the graveyard; the escape is
undecided, exactly Part III's 'live possibility of being wrong'.
"""

from __future__ import annotations

import numpy as np

S_LEP_BOUND = 0.1


def spectral_V(s, M, Nc=3):
    """Vector spectral function of a fermion of mass M (zero below threshold)."""
    s = np.asarray(s, dtype=float)
    out = np.zeros_like(s)
    m = s > 4 * M**2
    x = 4 * M**2 / s[m]
    out[m] = (Nc / (12 * np.pi)) * (1 + 0.5 * x) * np.sqrt(1 - x)
    return out


def spectral_A(s, M, Nc=3):
    """Axial-vector spectral function of a fermion of mass M (transverse continuum)."""
    s = np.asarray(s, dtype=float)
    out = np.zeros_like(s)
    m = s > 4 * M**2
    x = 4 * M**2 / s[m]
    out[m] = (Nc / (12 * np.pi)) * (1 - x) * np.sqrt(1 - x)
    return out


def spectral_difference(s, M, Nc=3):
    """rho_V - rho_A = (Nc/12pi)(6 M^2/s) sqrt(1-4M^2/s): the chiral-odd, positive,
    UV-finite difference that drives S > 0."""
    return spectral_V(s, M, Nc) - spectral_A(s, M, Nc)


def s_leading(Nc=3, n_doublets=1.0):
    """Leading-order (constituent-loop) S of chirally-broken doublets:
    S = N_c n_doublets / (6 pi)  (Peskin-Takeuchi)."""
    return Nc * n_doublets / (6.0 * np.pi)


def s_walking(fpi_over_MV, MA_over_MV=1.4):
    """The resonance-saturated (Weinberg-sum-rule) S that walking must achieve:
    S = 4 pi (f_pi/M_V)^2 (1 + M_V^2/M_A^2). S < 0.1 needs f_pi/M_V <~ 0.075."""
    x = 1.0 / MA_over_MV**2
    return 4.0 * np.pi * fpi_over_MV**2 * (1.0 + x)


def verdict(Nc=3, n_doublets=1.0):
    """Leading S vs the LEP bound, and the walking factor needed to escape."""
    S0 = s_leading(Nc, n_doublets)
    # f_pi/M_V needed for S < 0.1 (at QCD-like M_A/M_V)
    need = np.sqrt(S_LEP_BOUND / (4 * np.pi * (1 + 1 / 1.4**2)))
    return {
        "S_leading": float(S0),
        "in_graveyard": bool(S0 > S_LEP_BOUND),
        "fpi_over_MV_needed": float(need),
        "MV_over_fpi_needed": float(1 / need),
        "qcd_MV_over_fpi": 776.0 / 93.0,           # ~8.3
    }
