r"""A finite Woods-Saxon "bag" model of the torsiton tower -- the shape, not the exact solution.

The torsiton's effective potential is a FINITE bag well (Chapter 11; confirmed by the lattice that
the baryon binds and by the mean-field that V_eff = M(r)^2 rises to a dissociation threshold). NOT
Coulomb (no 1/r, so not the infinite Rydberg tower), NOT Morse (no hard inner wall). The clean
analytic stand-in is the Woods-Saxon well

    V(r) = -V0 / (1 + exp((r - R0)/a)),

a flat chiral-restored floor of depth V0, radius R0, surface thickness a. A finite well binds a
FINITE number of s-wave rungs -- the generations -- with the textbook hydrogen-like node structure
(0, 1, 2, ... interior nodes), and the count is set by the depth x width (sqrt(2 m V0) R0 / pi ~ N):
the cap. At V0=5, R0=3, a=0.5 it binds exactly three; a fourth appears by V0~6. The mass of each rung
is its overlap with the LOCAL dynamical mass M(r) (Eq. configmass weights the scalar density by M, not
by the condensate sigma) -- and M is large in the VACUUM, small in the chiral-restored core. So the
mass RISES with generation: the ground rung is tucked in the core (lightest), higher rungs spread into
the vacuum (heavier). The figure shows the wavefunctions, sigma(r), M(r), and the overlap-per-rung.

This is a MODEL (a single fermion in a fixed well), deliberately: it makes the shape, the node
structure, the cap, and the overlap mechanism legible. The depth and width are what the lattice
pins (depth ~ constituent mass^2 from m_N, width ~ the bag size from the nucleon's extent)."""
from __future__ import annotations

import numpy as np
from scipy.linalg import eigh_tridiagonal

_trapz = getattr(np, "trapezoid", None) or np.trapz


def ws_potential(r, V0, R0, a):
    """Woods-Saxon well V(r) = -V0/(1+exp((r-R0)/a)): floor -V0, radius R0, surface thickness a."""
    return -V0 / (1.0 + np.exp((np.asarray(r, dtype=float) - R0) / a))


def ws_spectrum(V0=5.0, R0=3.0, a=0.5, rmax=16.0, N=1600, m=1.0):
    """Bound s-wave levels of -(1/2m) u'' + V(r) u = E u in the Woods-Saxon well (reduced radial u,
    u(0)=u(rmax)=0). Returns (r, [(E_n, u_n), ...]) for the bound states E<0, normalised to
    int u_n^2 dr = 1, ground first. The number of them is the rung count (the generation cap)."""
    r = np.linspace(rmax / N, rmax, N)
    h = r[1] - r[0]
    V = ws_potential(r, V0, R0, a)
    diag = 1.0 / (m * h**2) + V
    off = -1.0 / (2.0 * m * h**2) * np.ones(N - 1)
    w, v = eigh_tridiagonal(diag, off)
    out = []
    for k in range(N):
        if w[k] < 0.0:
            u = v[:, k] / np.sqrt(h)
            if u[np.argmax(np.abs(u))] < 0:
                u = -u                                    # fix sign: positive first lobe
            out.append((float(w[k]), u))
    return r, out


def interior_nodes(u):
    """Number of interior sign changes of the reduced wavefunction u (0 for the ground rung)."""
    s = np.sign(u[np.abs(u) > 1e-6 * np.max(np.abs(u))])
    return int(np.sum(np.diff(s) != 0))


def overlap_masses(r, states, weight):
    """Configurational mass of each rung as the density overlap with a radial weight:
    m_n = int u_n(r)^2 weight(r) dr. For the physical (rising) ordering pass the LOCAL dynamical mass
    M(r) = 1 - sigma(r) (Eq. configmass: large in the vacuum); passing sigma(r) gives the wrong
    (falling) order. Returns the array (same order as states)."""
    w = np.asarray(weight, dtype=float)
    return np.array([_trapz(u**2 * w, r) for _, u in states])
