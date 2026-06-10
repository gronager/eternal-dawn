r"""Species-resolved genesis (Horizon, step 2) on a GPU-ready, BATCHED Skyrme engine.

The single-field genesis gives one species (baryons). Here several chiral SECTORS -- the fermion
towers, distinguished by their charge/colour coupling -- cool through ONE shared transition and each
condenses its OWN solitons, with a mass set by its coupling. A heavier coupling (a stiffer Skyrme
sector) makes a heavier knot, so the species fan out into a relative-mass hierarchy AS THEY CONDENSE
-- the dynamical version of horizon.py, where the masses come OUT of the soliton energies instead of
being read in.

THE ENGINE. All sectors are SU(2)/O(4) Skyrme fields evolved as ONE batched array n[B, L, L, L, 4]:
the batch axis B holds the species (and, if you like, seeds and cooling rates -- statistics and the
Kibble--Zurek scan for free). Every force/energy/baryon op is vectorised over B and over space, so the
whole ensemble advances in a single GPU pass. The backend is numpy (CPU) or cupy (GPU) by a one-line
swap -- the array API is identical; set CARTASIS_GPU=0 to force CPU. Time-stepping is sequential (a
Langevin Markov chain), but the per-step work -- all of it -- is parallel.

HONEST: a SCHEMATIC multi-sector model -- each tower an independent O(4) Skyrme sector with its
coupling from c_T -- not yet the full SU(3)xSU(2)xU(1) chiral field (that, with confinement binding
coloured quark-knots into colour-singlet baryons, is the next step). It shows the genuine new physics:
DISTINCT species condensing as distinct knots from one transition, masses ordered by their couplings.
"""
from __future__ import annotations

import os

# --- backend: cupy (GPU) if available and not disabled, else numpy (CPU). Identical array API. ---
if os.environ.get("CARTASIS_GPU", "1") != "0":
    try:
        import cupy as xp
        GPU = True
    except Exception:
        import numpy as xp
        GPU = False
else:
    import numpy as xp
    GPU = False

import numpy as np

_SP = (1, 2, 3)            # spatial axes of the batched field n[B, x, y, z, comp]


def _to_np(a):
    return xp.asnumpy(a) if GPU else np.asarray(a)


def _unit(n):
    return n / xp.sqrt(xp.sum(n * n, axis=-1, keepdims=True)).clip(1e-12)


def _d(f, i):     # central difference along spatial axis i (periodic)
    return 0.5 * (xp.roll(f, -1, i) - xp.roll(f, 1, i))


def _lap(f):      # 7-point periodic Laplacian over the 3 spatial axes
    return sum(xp.roll(f, s, ax) for ax in _SP for s in (-1, 1)) - 6.0 * f


def _cross4(u, v, w):
    """4D generalised cross product (u x v x w)_a = eps_{abcd} u_b v_c w_d, batched over leading axes."""
    out = xp.empty_like(u)
    for a in range(4):
        keep = [c for c in range(4) if c != a]
        p, q, r = keep
        det = (u[..., p] * (v[..., q] * w[..., r] - v[..., r] * w[..., q])
               - u[..., q] * (v[..., p] * w[..., r] - v[..., r] * w[..., p])
               + u[..., r] * (v[..., p] * w[..., q] - v[..., q] * w[..., p]))
        out[..., a] = (1.0 if a % 2 == 0 else -1.0) * det
    return out


_BNORM = 1.0 / (2.0 * np.pi ** 2)


def baryon_density(n):
    """Per-sector baryon density b[B, x, y, z] (x _BNORM integrates to the integer baryon number)."""
    A, Bv, C = _d(n, 1), _d(n, 2), _d(n, 3)
    return xp.sum(n * _cross4(A, Bv, C), axis=-1) * _BNORM


def _bcast(c):    # per-batch coupling (B,) -> (B,1,1,1,1) for broadcasting against n
    return xp.asarray(c, dtype=float).reshape((-1, 1, 1, 1, 1))


def energy_per_sector(n, kappa, mu2):
    """Total Skyrme energy of each sector: returns a length-B array (sigma + Skyrme + potential)."""
    D = [_d(n, i) for i in _SP]
    e = 0.5 * sum(xp.sum(d * d, -1) for d in D)
    sk = 0.0
    for a in range(3):
        for b in range(a + 1, 3):
            sk = sk + (xp.sum(D[a] * D[a], -1) * xp.sum(D[b] * D[b], -1)
                       - xp.sum(D[a] * D[b], -1) ** 2)
    kap = xp.asarray(kappa, dtype=float).reshape((-1, 1, 1, 1))
    mu = xp.asarray(mu2, dtype=float).reshape((-1, 1, 1, 1))
    e = e + 0.5 * kap * sk + mu * (1.0 - n[..., 0])
    return _to_np(xp.sum(e, axis=_SP))


def force(n, kappa, mu2):
    """Tangential force on each sector (batched): sigma + Skyrme W-current + potential, projected."""
    D = [_d(n, i) for i in _SP]
    F = _lap(n)
    kap = _bcast(kappa)
    for ii, i in enumerate(_SP):
        Wi = xp.zeros_like(n)
        for jj, j in enumerate(_SP):
            if j == i:
                continue
            dj2 = xp.sum(D[jj] * D[jj], -1, keepdims=True)
            dij = xp.sum(D[ii] * D[jj], -1, keepdims=True)
            Wi = Wi + dj2 * D[ii] - dij * D[jj]
        F = F + kap * _d(Wi, i)
    F[..., 0] += _bcast(mu2)[..., 0]
    return F - xp.sum(F * n, axis=-1, keepdims=True) * n


def run_species(species, L=40, steps=1500, dt=0.006, gamma=0.5, T0=0.8, T1=0.0, seed=0, record=8):
    """Cool B chiral SECTORS through one shared transition and watch each condense its species.

    `species` is a list of (name, kappa, mu2). All sectors are batched into n[B, L, L, L, 4] and
    advanced together (one GPU pass per step). Returns per-sector time series: energy, baryon content
    sum|b|, net B, the temperature schedule, the final baryon-density volumes, and the relic energy
    (the species' relative MASS). dt must suit the STIFFEST (largest kappa) sector."""
    names = [s[0] for s in species]
    kappa = np.array([s[1] for s in species], dtype=float)
    mu2 = np.array([s[2] for s in species], dtype=float)
    B = len(species)
    rng = xp.random.default_rng(seed) if GPU else np.random.default_rng(seed)
    n = _unit(rng.standard_normal((B, L, L, L, 4)))
    pi = xp.zeros_like(n)

    energy, content, netB, temps = [], [], [], []
    rec_at = set(np.linspace(0, steps - 1, record).astype(int))
    frames = None
    for s in range(steps):
        T = T0 + (T1 - T0) * (s / max(steps - 1, 1))
        F = force(n, kappa, mu2)
        xi = xp.sqrt(2.0 * gamma * max(T, 0.0) * dt) * rng.standard_normal(n.shape)
        pi = pi + dt * F - dt * gamma * pi + xi
        pi = pi - xp.sum(pi * n, axis=-1, keepdims=True) * n
        n = _unit(n + dt * pi)
        if s in rec_at:
            b = baryon_density(n)
            energy.append(energy_per_sector(n, kappa, mu2))
            content.append(_to_np(xp.sum(xp.abs(b), axis=_SP)))
            netB.append(_to_np(xp.sum(b, axis=_SP)))
            temps.append(T)
            frames = _to_np(b)
    energy = np.array(energy)                              # (record, B)
    return {"names": names, "kappa": kappa, "mu2": mu2, "energy": energy,
            "content": np.array(content), "netB": np.array(netB), "T": np.array(temps),
            "final_b": frames, "gpu": GPU, "L": L}


def _hedgehog_batch(B, L, width=5.0):
    """B copies of a B=1 hedgehog (the sectors differ by coupling, not by initial condition)."""
    g = [np.arange(L) - L / 2.0 for _ in range(3)]
    X, Y, Z = np.meshgrid(*g, indexing="ij")
    r = np.sqrt(X * X + Y * Y + Z * Z).clip(1e-6)
    F = np.pi * np.exp(-r / width)
    sF = np.sin(F)
    n1 = np.stack([np.cos(F), sF * X / r, sF * Y / r, sF * Z / r], axis=-1)
    return xp.asarray(np.broadcast_to(n1, (B,) + n1.shape).copy())


def relax_solitons(species, L=40, steps=500, dt=0.005, gamma=0.8, width=5.0):
    """The clean species MASS: relax ONE B=1 knot per sector (batched) at zero temperature and read
    its energy. `species` = [(name, kappa, mu2), ...]. Returns per-sector mass (energy), the relative
    mass (lightest=1), and the baryon number (should stay ~1). This is the species rest mass -- the
    masses come OUT of the soliton dynamics given the couplings (M ~ sqrt(kappa) for the Skyrme knot),
    the dynamical core of the Horizon spectrum."""
    names = [s[0] for s in species]
    kappa = np.array([s[1] for s in species], dtype=float)
    mu2 = np.array([s[2] for s in species], dtype=float)
    n = _unit(_hedgehog_batch(len(species), L, width))
    pi = xp.zeros_like(n)
    for _ in range(steps):
        F = force(n, kappa, mu2)
        pi = pi + dt * F - dt * gamma * pi
        pi = pi - xp.sum(pi * n, axis=-1, keepdims=True) * n
        n = _unit(n + dt * pi)
    mass = energy_per_sector(n, kappa, mu2)
    Bnum = _to_np(xp.sum(baryon_density(n), axis=_SP))
    return {"names": names, "kappa": kappa, "mass": mass, "mass_rel": mass / mass.min(),
            "B": Bnum, "gpu": GPU}
