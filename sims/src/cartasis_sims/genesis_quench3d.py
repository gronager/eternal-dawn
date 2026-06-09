r"""Genesis in 3D: torsitons (baryons) condensing from the cooling chiral fluid -- the real Skyrme model.

The 2D baby-Skyrme prototype proved the mechanism; this is the quantitative engine. The chiral
condensate is now the SU(2) field U = n0 + i n.tau, i.e. a 4-component unit vector n in S^3 on a 3D
grid -- the *standard nuclear Skyrme model*, which gives baryon masses and radii to ~30% with the
chiral scale as its only input. The energy

    E = int d^3x [ 1/2 sum_i |d_i n|^2  +  kappa/2 sum_{i<j} (|d_i n|^2|d_j n|^2 - (d_i n . d_j n)^2)
                   +  mu2 (1 - n0) ],

is the sigma gradient + the Skyrme stabiliser (finite knot size, Derrick) + a vacuum potential. The
conserved BARYON number is the degree of the map n: R^3 -> S^3,

    B = c int d^3x  eps^{ijk} eps_{abcd} n_a (d_i n_b)(d_j n_c)(d_k n_d),

normalised so a single hedgehog skyrmion gives B = 1. So in 3D the topological knots are literally
baryons (the torsitons), and a cooling quench condenses them out of the fluid -- mass = the knot's
field energy (the dominant ~94% of the configurational mass), abundance = the Kibble--Zurek defect
density (how much matter a cooling bounce makes). Pin kappa, mu2 to the lattice scale and the knot
energy becomes a torsiton mass in physical units (order of magnitude).

HONEST: still classical stochastic EFT (formation, mass to ~30%, abundance), not the exact quantum
spectrum (the lattice). For production (128^3) swap numpy -> cupy/jax: the array API is identical.
This module validates the physics on modest grids."""
from __future__ import annotations

import numpy as np

_DIRS = (0, 1, 2)


def _unit(n):
    return n / np.sqrt(np.sum(n * n, axis=-1, keepdims=True)).clip(1e-12)


def _d(f, i):  # central difference along spatial axis i (periodic)
    return 0.5 * (np.roll(f, -1, i) - np.roll(f, 1, i))


def _lap(f):  # 7-point periodic Laplacian in 3D
    return (sum(np.roll(f, s, ax) for ax in _DIRS for s in (-1, 1)) - 6.0 * f)


def _cross4(u, v, w):
    """The 4D generalised cross product (u x v x w)_a = eps_{abcd} u_b v_c w_d: the 4-vector orthogonal
    to u,v,w, via cofactor 3x3 dets of the (u,v,w) columns with axis a removed (signs (+,-,+,-))."""
    out = np.empty_like(u)
    cols = [0, 1, 2, 3]
    for a in range(4):
        keep = [c for c in cols if c != a]
        p, q, r = keep
        det = (u[..., p] * (v[..., q] * w[..., r] - v[..., r] * w[..., q])
               - u[..., q] * (v[..., p] * w[..., r] - v[..., r] * w[..., p])
               + u[..., r] * (v[..., p] * w[..., q] - v[..., q] * w[..., p]))
        out[..., a] = (1.0 if a % 2 == 0 else -1.0) * det
    return out


def baryon_density(n):
    """Unnormalised baryon (topological) density n . (d_x n x d_y n x d_z n); sum x B-norm = B in Z."""
    A, B, C = _d(n, 0), _d(n, 1), _d(n, 2)
    return np.sum(n * _cross4(A, B, C), axis=-1)


# normalisation so a single hedgehog skyrmion integrates to B = 1 (the 1/(2 pi^2) of the degree map,
# pinned numerically on a reference hedgehog -- see hedgehog_skyrmion / total_baryon).
_BNORM = 1.0 / (2.0 * np.pi ** 2)


def total_baryon(n):
    """The integer baryon number B = sum_x baryon_density * norm -- the conserved torsiton tally."""
    return float(baryon_density(n).sum() * _BNORM)


def energy_density(n, kappa=30.0, mu2=0.008):
    """Skyrme energy density: sigma gradient + Skyrme stabiliser (3 spatial pairs) + vacuum potential."""
    D = [_d(n, i) for i in _DIRS]
    e_sigma = 0.5 * sum(np.sum(d * d, -1) for d in D)
    e_sk = 0.0
    for i in _DIRS:
        for j in _DIRS:
            if j <= i:
                continue
            di2 = np.sum(D[i] * D[i], -1)
            dj2 = np.sum(D[j] * D[j], -1)
            dij = np.sum(D[i] * D[j], -1)
            e_sk += di2 * dj2 - dij * dij
    e_sk *= 0.5 * kappa
    e_pot = mu2 * (1.0 - n[..., 0])
    return e_sigma + e_sk + e_pot


def force(n, kappa=30.0, mu2=0.008):
    """Tangential force F = -dE/dn on S^3 (sigma + Skyrme W-current + potential), projected off n."""
    D = [_d(n, i) for i in _DIRS]
    F = _lap(n)                                              # sigma: +lap n
    for i in _DIRS:                                          # Skyrme: +kappa sum_i d_i W_i
        Wi = np.zeros_like(n)
        for j in _DIRS:
            if j == i:
                continue
            dj2 = np.sum(D[j] * D[j], -1, keepdims=True)
            dij = np.sum(D[i] * D[j], -1, keepdims=True)
            Wi += dj2 * D[i] - dij * D[j]
        F += kappa * _d(Wi, i)
    F[..., 0] += mu2                                         # potential pulls n0 -> +1 (vacuum)
    return F - np.sum(F * n, axis=-1, keepdims=True) * n     # tangent projection on S^3


def hedgehog_skyrmion(L, width=4.0, center=None):
    """A B=1 hedgehog: U = exp(i F(r) rhat.tau), n = (cos F, sin F rhat), F(0)=pi, F(inf)=0."""
    if center is None:
        center = (L / 2.0,) * 3
    g = [np.arange(L) - center[k] for k in range(3)]
    X, Y, Z = np.meshgrid(g[0], g[1], g[2], indexing="ij")
    r = np.sqrt(X * X + Y * Y + Z * Z).clip(1e-6)
    F = np.pi * np.exp(-r / width)
    s = np.sin(F)
    n = np.stack([np.cos(F), s * X / r, s * Y / r, s * Z / r], axis=-1)
    return _unit(n)


def run_quench3d(L=48, steps=2000, dt=0.006, kappa=30.0, mu2=0.008, gamma=0.5,
                 T0=0.8, T1=0.0, seed=0, record=10, hot_start=True, n_init=None):
    """Cool a hot, disordered SU(2) chiral fluid through the transition and watch BARYONS (torsitons)
    condense (underdamped Langevin, linear T ramp). Returns recorded baryon-density volumes, the net
    B, the unsigned content sum|B|, the energy, and the temperature -- the 3D genesis."""
    rng = np.random.default_rng(seed)
    if n_init is not None:
        n = _unit(np.array(n_init, dtype=float))
    elif hot_start:
        n = _unit(rng.standard_normal((L, L, L, 4)))
    else:
        n = _unit(np.tile([1.0, 0, 0, 0], (L, L, L, 1)) + 0.3 * rng.standard_normal((L, L, L, 4)))
    pi = np.zeros_like(n)
    frames, Bs, content, energy, temps = [], [], [], [], []
    rec_at = set(np.linspace(0, steps - 1, record).astype(int))
    for s in range(steps):
        T = T0 + (T1 - T0) * (s / max(steps - 1, 1))
        F = force(n, kappa, mu2)
        xi = np.sqrt(2.0 * gamma * max(T, 0.0) * dt) * rng.standard_normal(n.shape)
        pi = pi + dt * F - dt * gamma * pi + xi
        pi = pi - np.sum(pi * n, axis=-1, keepdims=True) * n
        n = _unit(n + dt * pi)
        if s in rec_at:
            b = baryon_density(n) * _BNORM
            frames.append((s, b.copy()))
            Bs.append(float(b.sum()))
            content.append(float(np.abs(b).sum()))
            energy.append(float(energy_density(n, kappa, mu2).sum()))
            temps.append(T)
    return {"frames": frames, "B": np.array(Bs), "content": np.array(content),
            "energy": np.array(energy), "T": np.array(temps), "L": L,
            "params": dict(kappa=kappa, mu2=mu2, gamma=gamma, dt=dt, steps=steps)}
