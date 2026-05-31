r"""Axis-alignment analysis (Phase-0, Corner B).

Tests the Cartasis prediction that a *single* axis underlies the large-angle CMB
anomalies and the galaxy-rotation asymmetry -- the imprint of a rotating parent
black hole's spin. We work with *axes* (headless directions on the sphere,
RP^2): each anomaly defines a line, not an arrow.

Two statistics:

* pairwise acute angle between two axes, with the isotropic p-value
  p(<= theta) = 1 - cos(theta)  (probability two random axes fall within the
  acute angle theta of each other);
* a concentration test for N axes: the largest eigenvalue tau1 of the
  orientation tensor T = (1/N) sum_i x_i x_i^T, calibrated by Monte Carlo
  against isotropically distributed axes.

EPISTEMIC HEALTH WARNING. These axes were noticed *because* they look aligned,
so naive p-values are a-posteriori-inflated lower bounds. Worse, the galaxy-spin
axis reportedly tracks the *Galactic* poles -- exactly what a Milky-Way
observational systematic (dust, inclination-classification bias) would produce.
The cosmological (parent-spin) interpretation only survives if the galaxy-spin
axis aligns with the *CMB* axis better than it aligns with the Galactic pole.
This module is built to make that comparison explicit, not to manufacture
significance.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


# ---------------------------------------------------------------------------
# Curated axis directions (Galactic coordinates, degrees).
# Values flagged verify=True still need checking against the primary source.
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Axis:
    name: str
    l: float            # Galactic longitude, deg
    b: float            # Galactic latitude, deg
    kind: str           # 'cmb' | 'galaxy' | 'geometry'
    reference: str
    verify: bool = True
    note: str = ""


AXES: list[Axis] = [
    Axis("CMB axis of evil (l=2/3)", 260.0, 60.0, "cmb",
         "Land & Magueijo 2005; Copi et al. review", verify=True,
         note="combined quadrupole-octupole preferred axis, ~(260,60)"),
    Axis("CMB quadrupole (l=2)", 240.0, 63.0, "cmb",
         "de Oliveira-Costa et al. 2004", verify=True),
    Axis("CMB octupole (l=3)", 308.0, 63.0, "cmb",
         "de Oliveira-Costa et al. 2004", verify=True),
    Axis("CMB kinematic dipole", 264.02, 48.25, "geometry",
         "Planck 2018", verify=False,
         note="our peculiar motion; partly a frame effect, not an anomaly"),
    Axis("CMB cold spot", 209.0, -57.0, "cmb",
         "Cruz et al. 2005", verify=True),
    Axis("North ecliptic pole", 96.38, 29.81, "geometry",
         "fixed Solar-System geometry", verify=False),
    Axis("Galactic pole", 0.0, 90.0, "geometry",
         "definition", verify=False,
         note="Milky-Way axis; the chief systematic confound"),
    Axis("Galaxy-spin asymmetry", 0.0, 90.0, "galaxy",
         "Shamir 2025 (JADES); DESI all-sky peaks at Galactic poles",
         verify=True,
         note="reported to peak toward the Galactic poles -- "
              "encoded as the Galactic-pole axis pending an all-sky dipole fit"),
]


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
def lb_to_vec(l_deg: float, b_deg: float) -> np.ndarray:
    """Galactic (l, b) in degrees -> unit vector in Galactic Cartesian frame."""
    l = np.radians(l_deg)
    b = np.radians(b_deg)
    return np.array([np.cos(b) * np.cos(l),
                     np.cos(b) * np.sin(l),
                     np.sin(b)])


def acute_angle_deg(v1: np.ndarray, v2: np.ndarray) -> float:
    """Acute angle (0..90 deg) between two axes (sign-independent)."""
    c = np.clip(abs(np.dot(v1, v2)), 0.0, 1.0)
    return float(np.degrees(np.arccos(c)))


def axis_pvalue(theta_deg: float) -> float:
    """Isotropic probability that two random axes fall within acute angle theta.

    For headless axes, p(<= theta) = 1 - cos(theta), theta in [0, 90].
    """
    return float(1.0 - np.cos(np.radians(theta_deg)))


# ---------------------------------------------------------------------------
# Concentration of N axes
# ---------------------------------------------------------------------------
def orientation_tensor(vectors: np.ndarray) -> np.ndarray:
    """T = (1/N) sum x_i x_i^T  (3x3, sign-invariant -> proper for axes)."""
    V = np.atleast_2d(vectors)
    return (V.T @ V) / V.shape[0]


def concentration(vectors: np.ndarray):
    """Return (tau1, principal_axis): largest eigenvalue and its eigenvector.

    tau1 in [1/3, 1]; tau1 -> 1 means all axes share one direction, tau1 -> 1/3
    is isotropic.
    """
    T = orientation_tensor(vectors)
    w, Q = np.linalg.eigh(T)
    return float(w[-1]), Q[:, -1]


def _random_axes(n: int, rng: np.random.Generator) -> np.ndarray:
    g = rng.standard_normal((n, 3))
    return g / np.linalg.norm(g, axis=1, keepdims=True)


def monte_carlo_pvalue(vectors: np.ndarray, n_trials: int = 200_000,
                       seed: int = 0) -> dict:
    """p that N isotropic axes are at least as concentrated (tau1) as observed.

    Fully vectorised: draws (n_trials, N) isotropic axes, forms each trial's
    orientation tensor, and compares the largest eigenvalue to the observed one.
    """
    V = np.atleast_2d(vectors)
    n = V.shape[0]
    tau1_obs, axis = concentration(V)
    rng = np.random.default_rng(seed)
    g = rng.standard_normal((n_trials, n, 3))
    g /= np.linalg.norm(g, axis=2, keepdims=True)
    T = np.einsum("tij,tik->tjk", g, g) / n          # (n_trials, 3, 3)
    tau1 = np.linalg.eigvalsh(T)[:, -1]               # largest eigenvalue
    return {"n_axes": n, "tau1": tau1_obs,
            "p_value": float(np.mean(tau1 >= tau1_obs)),
            "principal_axis_vec": axis}
