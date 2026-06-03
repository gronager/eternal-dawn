r"""De-confounding the galaxy-spin axis (Corner B, follow-up).

Fits the dipole axis of a galaxy-rotation-direction catalogue ourselves and asks
the discriminating question: is the asymmetry a genuine sky dipole, and does its
axis favour the CMB axis or merely the Galactic pole (a Milky-Way systematic)?

Estimator. With spin s_i in {+1 (CW), -1 (CCW)} at sky unit vector n_i, the
dipole vector is D = <s_i n_i>. Its direction is the spin axis; for full sky the
amplitude is A ~ 3|D|. Significance of a dipole *beyond* the overall asymmetry
comes from permuting the spins (label shuffle) and from bootstrapping the axis.

Caveats are structural, not incidental: a partial-sky footprint biases D toward
the survey's geometry, and human perception bias (Land et al. 2008) can fake an
overall asymmetry. The functions here expose those, they do not hide them.
"""

from __future__ import annotations

import numpy as np

# ICRS/equatorial -> Galactic rotation matrix (Hipparcos/ICRS conventions).
_A_G = np.array([
    [-0.054875560416, -0.873437090235, -0.483835015548],
    [ 0.494109427875, -0.444829629960,  0.746982244497],
    [-0.867666149019, -0.198076373431,  0.455983776175],
])


def equatorial_to_galactic(ra_deg, dec_deg):
    """(RA, Dec) J2000 in degrees -> (l, b) Galactic in degrees."""
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)
    e = np.array([np.cos(dec) * np.cos(ra),
                  np.cos(dec) * np.sin(ra),
                  np.sin(dec)])
    g = _A_G @ e
    l = np.degrees(np.arctan2(g[1], g[0])) % 360.0
    b = np.degrees(np.arcsin(np.clip(g[2], -1.0, 1.0)))
    return l, b


def lb_to_vec(l_deg, b_deg):
    l = np.radians(l_deg)
    b = np.radians(b_deg)
    return np.array([np.cos(b) * np.cos(l),
                     np.cos(b) * np.sin(l),
                     np.sin(b)])


def _acute_deg(v1, v2):
    return float(np.degrees(np.arccos(np.clip(abs(np.dot(v1, v2)), 0.0, 1.0))))


def fit_dipole(vecs: np.ndarray, spin: np.ndarray) -> dict:
    """Dipole vector D = <s_i n_i>; returns its (l,b) direction and amplitude.

    `vecs` is (N,3) Galactic unit vectors, `spin` is (N,) in {+1,-1}.
    """
    D = (spin[:, None] * vecs).mean(axis=0)
    mag = float(np.linalg.norm(D))
    axis = D / mag if mag > 0 else np.array([0.0, 0.0, 1.0])
    l = np.degrees(np.arctan2(axis[1], axis[0])) % 360.0
    b = np.degrees(np.arcsin(np.clip(axis[2], -1.0, 1.0)))
    return {"D": D, "mag": mag, "amplitude_3D": 3.0 * mag,
            "axis": axis, "l": float(l), "b": float(b)}


def dipole_null_pvalue(vecs, spin, n_perm=2000, seed=0):
    """p that label-shuffled spins produce a dipole as large as observed."""
    rng = np.random.default_rng(seed)
    obs = np.linalg.norm((spin[:, None] * vecs).mean(axis=0))
    s = spin.copy()
    count = 0
    for _ in range(n_perm):
        rng.shuffle(s)
        count += np.linalg.norm((s[:, None] * vecs).mean(axis=0)) >= obs
    return float(count / n_perm)


def bootstrap_axis(vecs, spin, n_boot=2000, seed=1):
    """Bootstrap the dipole axis; return (mean_axis, 68%% angular radius deg)."""
    rng = np.random.default_rng(seed)
    n = len(spin)
    axes = np.empty((n_boot, 3))
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        D = (spin[idx, None] * vecs[idx]).mean(axis=0)
        axes[i] = D / np.linalg.norm(D)
    mean_axis = axes.mean(axis=0)
    mean_axis /= np.linalg.norm(mean_axis)
    ang = np.degrees(np.arccos(np.clip(np.abs(axes @ mean_axis), 0.0, 1.0)))
    return mean_axis, float(np.percentile(ang, 68.0))


def analyse(l, b, spin, cmb_axis_lb=(260.0, 60.0)):
    """Full de-confounding summary for a spin catalogue (Galactic l, b, spin)."""
    from scipy.stats import binomtest

    vecs = np.array([lb_to_vec(li, bi) for li, bi in zip(l, b)])
    spin = np.asarray(spin, dtype=float)
    n_cw = int((spin > 0).sum())
    n_ccw = int((spin < 0).sum())
    n = n_cw + n_ccw
    f_cw = n_cw / n
    binom_p = binomtest(n_cw, n, 0.5).pvalue

    dip = fit_dipole(vecs, spin)
    null_p = dipole_null_pvalue(vecs, spin)
    mean_axis, axis_68 = bootstrap_axis(vecs, spin)

    cmb = lb_to_vec(*cmb_axis_lb)
    gal_pole = np.array([0.0, 0.0, 1.0])
    ang_cmb = _acute_deg(dip["axis"], cmb)
    ang_pole = _acute_deg(dip["axis"], gal_pole)

    return {
        "n": n, "n_cw": n_cw, "n_ccw": n_ccw, "f_cw": f_cw,
        "excess_percent": 100.0 * (n_cw - n_ccw) / n, "binom_p": float(binom_p),
        "dipole_lb": (dip["l"], dip["b"]), "amplitude_3D": dip["amplitude_3D"],
        "dipole_null_p": null_p, "axis_68_deg": axis_68,
        "angle_to_cmb_axis": ang_cmb, "angle_to_galactic_pole": ang_pole,
        "vecs": vecs, "spin": spin,
    }


def load_spin_catalogue(path):
    """Load a galaxy-spin catalogue and return (l, b, spin) in Galactic degrees.

    Accepts two schemas:
      * GZ1-style   : columns RAJ2000, DEJ2000, pcS, paS (spin = +1 if pcS>paS);
      * generic     : columns ra_deg, dec_deg, spin (spin in {+1 CW, -1 CCW}).
    """
    d = np.genfromtxt(path, delimiter=",", names=True)
    cols = d.dtype.names
    if "pcS" in cols and "paS" in cols:
        ra, dec = d["RAJ2000"], d["DEJ2000"]
        spin = np.where(d["pcS"] > d["paS"], 1.0, -1.0)
    elif "spin" in cols:
        ra = d["ra_deg"] if "ra_deg" in cols else d["RAJ2000"]
        dec = d["dec_deg"] if "dec_deg" in cols else d["DEJ2000"]
        spin = d["spin"].astype(float)
    else:
        raise ValueError(f"unrecognised spin-catalogue schema: {cols}")
    lb = np.array([equatorial_to_galactic(r, c) for r, c in zip(ra, dec)])
    return lb[:, 0], lb[:, 1], spin
