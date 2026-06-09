r"""Horizon, step 3a: CONFINEMENT -- why we see baryons and leptons, not free quarks.

Step 2 condensed independent species. The new physics here is the one that distinguishes them: the
strong-torsion sector CONFINES (the framework's L1/L2: a linear quark potential V(r)=sigma r, with
sigma = 2 pi v^2 -- one substrate sets both the mass and the string tension). So a COLOURED knot (a
quark) cannot exist free -- its colour-electric string costs energy proportional to its separation --
while COLOUR-SINGLET combinations (mesons q-qbar, baryons q-q-q) and COLOURLESS knots (leptons) are
finite-energy and free. That is exactly the observed pattern: baryons and leptons, never a bare quark.

This module makes that quantitative and ties it to the lattice number sigma:
  * the confining potential V(r) = sigma r (the measured string),
  * the species ENERGETICS: a lepton is free (E=M0), a single quark is CONFINED (E grows with the box
    -- no finite free-quark state), a meson binds at a finite size, a baryon is finite,
  * a small N-body binding demo: coloured charges pulled together by the string into colour-singlet
    hadrons while colourless leptons diffuse free.

HONEST: a U(1)-colour PROXY of confinement (the exact-in-2D linear potential), not literal SU(3) flux
tubes -- the mechanism (gauge confinement binds charge into neutral composites), made tractable. The
full SU(3) flavour Skyrme baryon octet and the electroweak lepton sector are scoped in HORIZON.md."""
from __future__ import annotations

import numpy as np


def string_potential(r, sigma, core=0.5):
    """The colour-singlet (quark-antiquark) potential: a linear confining string V=sigma r plus a
    short-range core that gives a hadron a finite SIZE. Minimised at r* = sqrt(core/sigma)."""
    r = np.asarray(r, dtype=float)
    return sigma * r + core / np.clip(r, 1e-9, None)


def meson_size(sigma, core=0.5):
    """The bound size r* of a q-qbar meson: minimise sigma r + core/r."""
    return float(np.sqrt(core / sigma))


def species_energy(kind, sigma, M0=1.0, box=20.0, core=0.5):
    """Rest + confinement energy of a species (in units of the knot mass M0). The point: a 'quark'
    has NO finite free-state energy (it grows with the box -- confined), while 'lepton', 'meson',
    'baryon' are finite and box-independent.
      kind in {'lepton', 'quark', 'meson', 'baryon'}."""
    if kind == "lepton":                                   # colourless: free
        return M0
    if kind == "quark":                                    # isolated colour: string to the boundary
        return M0 + sigma * box                            # diverges with the box -> CONFINED
    if kind == "meson":                                    # q-qbar: a string of the bound length
        rstar = meson_size(sigma, core)
        return 2.0 * M0 + sigma * rstar + core / rstar     # = 2 M0 + 2 sqrt(sigma core), finite
    if kind == "baryon":                                   # q-q-q colour singlet: a finite Y-string
        rstar = meson_size(sigma, core)
        return 3.0 * M0 + np.sqrt(3.0) * sigma * rstar     # finite (Y-junction ~ sqrt3 x leg)
    raise ValueError(f"unknown species {kind!r}")


def is_confined(kind, sigma, M0=1.0, core=0.5):
    """A species is CONFINED if its energy grows without bound as the box grows (no free state)."""
    e_small = species_energy(kind, sigma, M0=M0, box=10.0, core=core)
    e_big = species_energy(kind, sigma, M0=M0, box=1000.0, core=core)
    return bool(e_big > e_small + 1e-6)


def run_binding(n_quark_pairs=8, n_lepton=6, L=40.0, steps=4000, dt=0.02, sigma=0.08, core=4.0,
                gamma=0.5, T0=0.4, T1=0.0, seed=0):
    """A 2D N-body demo: coloured charges (+/-1, quarks/antiquarks) bound by the confining string into
    colour-singlet mesons, while colourless leptons (q=0) diffuse free. Underdamped Langevin with
    cooling. Returns the trajectory endpoints and per-species clustering, illustrating the pattern --
    bound hadrons, free leptons. (Mesons here for tractability; baryons = SU(3) Y-strings, HORIZON.md.)"""
    rng = np.random.default_rng(seed)
    q = np.concatenate([np.ones(n_quark_pairs), -np.ones(n_quark_pairs), np.zeros(n_lepton)])
    N = len(q)
    x = rng.uniform(0, L, (N, 2))
    v = np.zeros((N, 2))
    rstar = meson_size(sigma, core)

    def forces(x):
        d = x[:, None, :] - x[None, :, :]
        d -= L * np.round(d / L)                            # minimum image (periodic)
        r = np.sqrt((d * d).sum(-1)) + 1e-9
        rhat = d / r[..., None]
        qq = q[:, None] * q[None, :]
        # opposite colour -> confining string (attractive, linear) + core repulsion; same -> repel
        fmag = np.where(qq < 0, -sigma + core / r ** 2, np.where(qq > 0, core / r ** 2, 0.0))
        np.fill_diagonal(fmag, 0.0)
        fmag *= np.exp(-r / (6 * rstar))                   # local string (smooth cutoff, keeps it cheap)
        return (fmag[..., None] * rhat).sum(1)

    for s in range(steps):
        T = T0 + (T1 - T0) * (s / max(steps - 1, 1))
        F = forces(x)
        v = v + dt * F - dt * gamma * v + np.sqrt(2 * gamma * max(T, 0) * dt) * rng.standard_normal((N, 2))
        x = (x + dt * v) % L

    # nearest opposite-charge distance per coloured particle (bound mesons sit at ~r*)
    d = x[:, None, :] - x[None, :, :]
    d -= L * np.round(d / L)
    r = np.sqrt((d * d).sum(-1))
    np.fill_diagonal(r, np.inf)
    opp = (q[:, None] * q[None, :]) < 0
    nearest_opp = np.where(opp.any(1), np.where(opp, r, np.inf).min(1), np.nan)
    coloured = q != 0
    bound_frac = float(np.mean(nearest_opp[coloured] < 3 * rstar))
    return {"x": x, "q": q, "nearest_opp": nearest_opp, "r_star": rstar,
            "bound_fraction": bound_frac, "L": L, "sigma": sigma}
