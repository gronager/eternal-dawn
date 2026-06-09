r"""Genesis at the particle scale: torsitons condensing from a cooling chiral fluid (2D prototype).

THE PICTURE (the user's, made runnable). Just post-bounce the chiral order parameter is a hot,
disordered "tally" -- a relativistic quantum fluid with no preferred direction. As temperature drops
through the symmetry-breaking transition the fluid CONDENSES: domains order, and the conserved
topological winding gets trapped into localized spinning knots -- the torsitons. This is the
Kibble--Zurek mechanism, the same physics that freezes cosmic strings and superfluid vortices.

THE MODEL. The bosonized chiral condensate is an O(3) order parameter n(x) in S^2 (the "Fierzed
spin^2 term" becomes its potential). The energy is the baby-Skyrme functional

    E = int d^2x [ 1/2 |grad n|^2  +  kappa/2 |dx n x dy n|^2  +  mu2 (1 - n_z) ],

the sigma gradient (couples cells -- the "non-local" interaction), the Skyrme stabiliser (gives the
knot a FINITE size -- Derrick), and a potential pinning the vacuum to n=zhat. The conserved
"baryon/fermion number" is the topological charge

    Q = 1/4pi int d^2x  n . (dx n x dy n)  in Z,

the winding S^2->S^2 -- so the fluid's tally is carried TOPOLOGICALLY by the knots, conserved through
the quench. We evolve underdamped Langevin (the fluid has momentum AND dissipation -- a thermostat)
and ramp the temperature down: the field waves, then condenses, and the torsitons light up as peaks
of the charge density.

HONEST: this is CLASSICAL stochastic effective field theory -- it gives the FORMATION dynamics, the
defect content, and the Kibble--Zurek scaling, NOT the exact quantum masses (that is the lattice).
The parameters (kappa, mu2, the cooling) map to the microscopic torsion couplings -- run dimensionless
to see the genesis; pin the numbers from the lattice/mean-field later."""
from __future__ import annotations

import numpy as np

_TWO_PI4 = 1.0 / (4.0 * np.pi)


def _unit(n):
    """Project a 3-vector field (..., 3) back onto the unit sphere S^2."""
    return n / np.sqrt(np.sum(n * n, axis=-1, keepdims=True)).clip(1e-12)


def _dx(f):  # central difference along x (axis 0), periodic
    return 0.5 * (np.roll(f, -1, 0) - np.roll(f, 1, 0))


def _dy(f):  # central difference along y (axis 1), periodic
    return 0.5 * (np.roll(f, -1, 1) - np.roll(f, 1, 1))


def _lap(f):  # 5-point periodic Laplacian
    return (np.roll(f, -1, 0) + np.roll(f, 1, 0) + np.roll(f, -1, 1) + np.roll(f, 1, 1) - 4.0 * f)


def topological_charge_density(n):
    """q(x) = (1/4pi) n . (dx n x dy n) -- the local winding (torsiton) density; sums to Q in Z."""
    A, B = _dx(n), _dy(n)
    return _TWO_PI4 * np.sum(n * np.cross(A, B, axis=-1), axis=-1)


def total_charge(n):
    """The integer topological charge Q = sum_x q(x) -- the conserved fermion/baryon tally."""
    return float(topological_charge_density(n).sum())


def coarse_order(n, sigma=3.0):
    """The local condensate magnitude |<n>| (coarse-grained order parameter): 0 in the hot disordered
    fluid (the field cancels on averaging), ~1 in the cold condensed vacuum. This is the *medium* the
    torsitons are knotted in -- the chiral condensate itself, the thing the charge-density map omits."""
    from scipy.ndimage import gaussian_filter
    m = np.stack([gaussian_filter(n[..., a], sigma, mode="wrap") for a in range(n.shape[-1])], -1)
    return np.sqrt(np.sum(m * m, -1)).clip(0.0, 1.0)


def energy_density(n, kappa=0.5, mu2=0.1):
    """The baby-Skyrme energy density: sigma gradient + Skyrme stabiliser + vacuum potential."""
    A, B = _dx(n), _dy(n)
    e_sigma = 0.5 * (np.sum(A * A, -1) + np.sum(B * B, -1))
    c = np.cross(A, B, axis=-1)
    e_skyrme = 0.5 * kappa * np.sum(c * c, -1)
    e_pot = mu2 * (1.0 - n[..., 2])
    return e_sigma + e_skyrme + e_pot


def force(n, kappa=0.5, mu2=0.1):
    """The tangential force F = -dE/dn projected onto T_n S^2 (sigma + Skyrme + potential)."""
    A, B = _dx(n), _dy(n)
    c = np.cross(A, B, axis=-1)
    F = _lap(n)                                              # sigma:  +lap n
    F = F + kappa * (_dx(np.cross(B, c, axis=-1)) - _dy(np.cross(A, c, axis=-1)))   # Skyrme
    F[..., 2] += mu2                                         # potential pulls n_z -> +1 (vacuum)
    return F - np.sum(F * n, axis=-1, keepdims=True) * n     # project to the tangent plane


def baby_skyrmion(L, m=1, width=6.0, center=None):
    """A baby-Skyrmion ansatz of winding m: hedgehog profile f(0)=pi, f(inf)=0, charge Q = m."""
    if center is None:
        center = (L / 2.0, L / 2.0)
    x = np.arange(L) - center[0]
    y = np.arange(L) - center[1]
    X, Y = np.meshgrid(x, y, indexing="ij")
    r = np.sqrt(X * X + Y * Y)
    th = np.arctan2(Y, X)
    f = np.pi * np.exp(-r / width)                           # pi at core -> 0 outside
    n = np.stack([np.sin(f) * np.cos(m * th), np.sin(f) * np.sin(m * th), np.cos(f)], axis=-1)
    return _unit(n)


def run_quench(L=128, steps=6000, dt=0.08, kappa=0.6, mu2=0.12, gamma=0.6,
               T0=0.7, T1=0.0, seed=0, record=16, hot_start=True, n_init=None):
    """Cool a hot, disordered chiral fluid through the symmetry-breaking transition and watch the
    torsitons condense (underdamped Langevin, linear T ramp T0 -> T1).

    Returns a dict with the recorded snapshots (`frames`: list of (t, n, q)), the time series of the
    net charge `Q`, the unsigned winding content `content = sum|q|` (the torsiton count proxy, which
    falls as +/- pairs annihilate), the energy, and the temperature schedule. A faster ramp freezes in
    MORE torsitons -- the Kibble--Zurek signature (scan via `quench_rate_scan`)."""
    rng = np.random.default_rng(seed)
    if n_init is not None:
        n = _unit(np.array(n_init, dtype=float))
    elif hot_start:
        n = _unit(rng.standard_normal((L, L, 3)))            # disordered: the post-bounce "tally"
    else:
        n = _unit(np.tile([0.0, 0.0, 1.0], (L, L, 1)) + 0.3 * rng.standard_normal((L, L, 3)))
    pi = np.zeros_like(n)

    frames, Qs, content, energy, temps = [], [], [], [], []
    rec_at = set(np.linspace(0, steps - 1, record).astype(int))
    for s in range(steps):
        T = T0 + (T1 - T0) * (s / max(steps - 1, 1))         # linear cooling
        F = force(n, kappa, mu2)
        xi = np.sqrt(2.0 * gamma * max(T, 0.0) * dt) * rng.standard_normal(n.shape)
        pi = pi + dt * F - dt * gamma * pi + xi              # underdamped Langevin (momentum + noise)
        pi = pi - np.sum(pi * n, axis=-1, keepdims=True) * n  # keep momentum tangent
        n = _unit(n + dt * pi)                               # step + reproject to the sphere
        if s in rec_at:
            q = topological_charge_density(n)
            frames.append((s, n.copy(), q))
            Qs.append(float(q.sum()))
            content.append(float(np.abs(q).sum()))
            energy.append(float(energy_density(n, kappa, mu2).sum()))
            temps.append(T)
    return {"frames": frames, "Q": np.array(Qs), "content": np.array(content),
            "energy": np.array(energy), "T": np.array(temps), "L": L,
            "params": dict(kappa=kappa, mu2=mu2, gamma=gamma, dt=dt, steps=steps)}


def quench_rate_scan(rates, L=96, base_steps=8000, **kw):
    """Kibble--Zurek: final torsiton content vs cooling rate. `rates` multiply the cooling speed
    (fewer steps = faster quench). Returns (rates, final_content) -- faster quench freezes in more."""
    out = []
    for rate in rates:
        steps = max(400, int(base_steps / rate))
        res = run_quench(L=L, steps=steps, record=2, **kw)
        out.append(res["content"][-1])
    return np.array(rates, dtype=float), np.array(out)
