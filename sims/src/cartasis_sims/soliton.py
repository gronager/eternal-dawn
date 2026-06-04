r"""Part II/IV (gate calculation): the gravity--torsion soliton spectrum.

THE OBJECT. Part II proposes that elementary fermions are stationary, regular,
spin-supported solitons of the Einstein--Cartan field equations. Integrating out
torsion leaves a Dirac field with a four-fermion (Hehl--Datta) self-interaction
~ (psi-bar gamma5 gamma^mu psi)^2 -- so the binding "potential" is NOT external: it
is the well the fermion's own axial-spin density digs for itself. The problem is a
NONLINEAR, self-consistent Dirac equation (NJL / Soler class), not a fixed 1-D well.

WHAT THIS MODULE DOES (a first approximation, honestly labelled). It solves the
radial Dirac equation in a *given* scalar confining well S(r) -- the effective mass
M(r) = m + S(r) the four-fermion interaction generates -- and reads off the bound
energy levels. This is the first, non-self-consistent pass: take the well as fixed,
get the tower; the self-consistent refinement (the well sourced by the levels it
binds, with Pauli filling) is the next step and is where the real difficulty lives.

DERRICK. Why a soliton can exist at all in 3-D, where Derrick's theorem forbids
static scalar lumps: the Dirac (first-order spinor) kinetic term scales differently
from a scalar's, and the binding here balances two terms with different size-scaling
(inward curvature/attraction vs the outward rho^2 torsion wall) -- the same
two-terms-balance that makes the cosmological bounce. Soler-type nonlinear Dirac
solitons are known to exist in 3-D for exactly this reason.

SCALE CAVEAT. The bare Einstein--Cartan four-fermion coupling is ~ G/c^4 (Planck
suppressed), so a soliton bound by it alone sits near the Planck mass, not the
electron mass. Whether the particle-scale binding density equals the cosmological
rho_C is the open "density question" (Part II target 5). Here we therefore study the
spectrum SHAPE in dimensionless units (the Regge slope, level ratios, spin-orbit),
which is what feeds the electroweak S-parameter -- the absolute scale is set
separately and left open.

Verified against SymPy (radial reduction) and a known limit (test suite).
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp


def _well(r, kind: str, depth: float, width: float):
    """Scalar confining well S(r) (adds to the mass): the self-generated four-fermion
    potential, in dimensionless units. `depth` and `width` set its scale/shape."""
    x = r / width
    if kind == "linear":          # relativistic-string / flux-tube confinement
        return depth * x
    if kind == "harmonic":        # Dirac-oscillator-like
        return depth * x**2
    if kind == "bounce":          # rho^2 torsion wall vs softer attraction (anharmonic)
        # steep outer wall, soft floor -- the qualitative bounce-well shape
        return depth * (x**2 + 0.5 * x**4)
    raise ValueError(f"unknown well: {kind}")


def _radial_rhs(r, y, E, m, kappa, kind, depth, width):
    # u = r*G, v = r*F: the REDUCED radial functions (regular: u ~ r^{|kappa|}).
    u, v = y
    S = _well(r, kind, depth, width)
    rr = max(r, 1e-9)
    du = -(kappa / rr) * u + (E + m + S) * v
    dv = (kappa / rr) * v - (E - m - S) * u
    return [du, dv]


def _turning_point(E, m, kind, depth, width):
    """Outer classical turning radius where m + S(r) = E (S(r_t) = E - m)."""
    s = (E - m) / depth
    if s <= 0:
        return 0.3 * width
    if kind == "linear":
        x2 = s
    elif kind == "harmonic":
        x2 = s
    else:  # bounce: 0.5 x^4 + x^2 - s = 0  ->  x^2 = -1 + sqrt(1+2s)
        x2 = -1.0 + np.sqrt(1.0 + 2.0 * s)
    return width * np.sqrt(max(x2, 0.0))


def wavefunction(E, m=1.0, kappa=-1, kind="bounce", depth=6.0, width=1.0,
                 margin=1.6, n_pts=1400, cap=1e4):
    """Integrate (G,F) from the origin to just past the turning point; truncate at
    the onset of the exponentially-growing (forbidden) branch to stay finite.
    Returns (r, G, F)."""
    r_stop = _turning_point(E, m, kind, depth, width) + margin * width
    r0 = 1e-3 * width
    r = np.linspace(r0, r_stop, n_pts)
    # regular start: u = r^{|kappa|}, v -> 0  (u = rG, v = rF)
    y0 = [r0 ** abs(kappa), 0.0]
    sol = solve_ivp(_radial_rhs, (r0, r_stop), y0, t_eval=r,
                    args=(E, m, kappa, kind, depth, width),
                    rtol=1e-8, atol=1e-10, method="DOP853")
    u, v = sol.y[0], sol.y[1]
    G, F = u / sol.t, v / sol.t         # back to G = u/r, F = v/r
    big = np.abs(G) > cap
    if big.any():                       # cut at the forbidden-region blow-up onset
        i = int(np.argmax(big))
        return sol.t[:i], G[:i], F[:i]
    return sol.t, G, F


def _node_count(E, m, kappa, kind, depth, width):
    """Number of nodes (sign changes) of the large component in the allowed region.
    Increments by one at each bound-state eigenvalue."""
    _, G, _ = wavefunction(E, m, kappa, kind, depth, width)
    if G.size < 3:
        return 0
    s = np.sign(G)
    s = s[s != 0]
    return int(np.sum(s[1:] != s[:-1]))


def energy_levels(n_levels: int = 6, m: float = 1.0, kappa: int = -1,
                  kind: str = "bounce", depth: float = 6.0, width: float = 1.0,
                  E_max: float = 30.0, n_scan: int = 200):
    """Bound-state energies E_n (E > 0 tower): a Dirac fermion in the scalar
    confining well -- the soliton's internal/relative ladder. Located by where the
    node count increments (robust against the forbidden-region blow-up)."""
    Es = np.linspace(m + 1e-2, E_max, n_scan)
    nodes = np.array([_node_count(E, m, kappa, kind, depth, width) for E in Es])
    levels = []
    for k in range(n_levels):
        idx = np.where(nodes == k + 1)[0]            # first E with k+1 nodes
        if idx.size == 0:
            break
        i = idx[0]
        a, b = Es[i - 1], Es[i]                       # bracket the increment
        for _ in range(40):
            c = 0.5 * (a + b)
            if _node_count(c, m, kappa, kind, depth, width) >= k + 1:
                b = c
            else:
                a = c
        levels.append(0.5 * (a + b))
    return np.array(levels)


def regge_slope(levels: np.ndarray, m: float = 1.0):
    """Fit m^2_n = E_n^2 vs level index n: a linear (Regge) fit returns the slope.
    A relativistic confining ladder gives E^2 ~ linear in n -- the hadron Regge form."""
    n = np.arange(len(levels))
    A = np.vstack([n, np.ones_like(n)]).T
    slope, intercept = np.linalg.lstsq(A, levels**2, rcond=None)[0]
    resid = levels**2 - (slope * n + intercept)
    rms = float(np.sqrt(np.mean(resid**2))) / float(np.mean(levels**2))
    return float(slope), float(intercept), rms


def effective_mass_profile(r, m: float = 1.0, kind: str = "bounce",
                           depth: float = 6.0, width: float = 1.0):
    """M(r) = m + S(r): the effective mass the fermion sees -- the self-generated
    four-fermion 'potential' to plot."""
    return m + _well(np.asarray(r, dtype=float), kind, depth, width)


# ---------------------------------------------------------------------------
# Electroweak S parameter (Part III make-or-break), Peskin--Takeuchi.
# ---------------------------------------------------------------------------
def s_parameter(fpi_over_MV: float, MA_over_MV: float = 1.585):
    r"""Peskin--Takeuchi S from a resonance-saturated, Weinberg-sum-rule sector:

        S = 4 pi (f_pi/M_V)^2 (1 + (M_V/M_A)^2)        (SymPy-verified).

    A QCD-like sector (f_pi/M_V ~ 0.12, M_A/M_V ~ 1.6) gives S ~ 0.25 -- the
    technicolor 'graveyard'. S < 0.1 demands f_pi/M_V <~ 0.075, i.e. M_V/f_pi >~ 13
    (vs QCD's ~8): heavier resonances relative to f_pi, the hallmark of 'walking'
    (near-conformal, anharmonic) dynamics."""
    x = 1.0 / MA_over_MV**2
    import math
    return 4.0 * math.pi * fpi_over_MV**2 * (1.0 + x)


QCD_FPI_OVER_MV = 93.0 / 776.0      # ~0.12
QCD_MA_OVER_MV = 1230.0 / 776.0     # ~1.585
S_LEP_BOUND = 0.1                   # rough electroweak-precision ceiling


def fpi_over_MV_for_S(S_target: float = S_LEP_BOUND, MA_over_MV: float = QCD_MA_OVER_MV):
    """The f_pi/M_V (walking-ness) needed to reach a target S."""
    import math
    x = 1.0 / MA_over_MV**2
    return math.sqrt(S_target / (4.0 * math.pi * (1.0 + x)))
