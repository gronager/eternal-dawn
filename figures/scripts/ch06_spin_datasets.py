#!/usr/bin/env python3
r"""The shared-axis test on TWO independent galaxy-spin catalogues (Ch.6).

The framework predicts one inherited axis behind the CMB 'axis of evil' and the
galaxy-spin handedness. Here we fit the spin-dipole axis OURSELVES on two genuinely
independent, real catalogues and ask whether they agree with each other and with the
CMB axis -- or merely with the survey/Galactic geometry (a systematic):

  * GZ1  (Lintott+ 2011): 30,126 crowd-classified CW/ACW spirals, northern SDSS.
  * Iye  (Iye+ 2019):       530 expert S/Z-winding spirals, broader sky.

Result (honest): the two methods DISAGREE on the overall asymmetry (GZ1 a -3.9% CW
excess at binomial p~1e-11; Iye +1.9%, p~0.7), NEITHER shows a significant sky dipole
(label-shuffle p~0.6-0.8), and BOTH best-fit axes sit within ~7 deg of the Galactic
pole -- ~30 deg from the CMB axis. So today's handedness data is systematics-limited:
method-dependent perception bias plus footprint geometry, not a clean cosmic axis.
This neither confirms nor refutes the prediction (a low-spin parent predicts a quiet
sky too); it fixes what a real test needs -- an all-sky, mirror-validated catalogue
(Euclid, Rubin).

Renders figures/pdf/spin_datasets.pdf, writes sims/output/spin_datasets.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import galaxy_spins as gs

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
DATA = os.path.join(ROOT, "data", "galaxy_spins")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

CMB_AXIS = (260.0, 60.0)        # 'axis of evil', Galactic (l, b)
ECL_POLE = (96.4, 29.8)         # north ecliptic pole in Galactic coords

CATALOGUES = [
    ("GZ1 (crowd CW/ACW, $N$=30126)", "gz1_raw_debiased.csv", "C0"),
    ("Iye 2019 (expert S/Z, $N$=530)", "iye2019.csv", "C1"),
]


def _aitoff_xy(l_deg, b_deg):
    """Galactic (l,b) deg -> aitoff axes coords (l wrapped to [-180,180])."""
    l = ((np.asarray(l_deg) + 180.0) % 360.0) - 180.0
    return np.radians(l), np.radians(b_deg)


def _small_circle(l_deg, b_deg, radius_deg, n=120):
    """Points (l,b) deg on the spherical small circle of angular radius `radius_deg`
    about the axis (l,b) -- handles the poles correctly."""
    axis = gs.lb_to_vec(l_deg, b_deg)
    # an orthonormal basis perpendicular to the axis
    ref = np.array([0.0, 0.0, 1.0]) if abs(axis[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    u = np.cross(axis, ref); u /= np.linalg.norm(u)
    v = np.cross(axis, u)
    r = np.radians(radius_deg)
    th = np.linspace(0, 2 * np.pi, n)
    pts = (np.cos(r) * axis[:, None]
           + np.sin(r) * (np.cos(th) * u[:, None] + np.sin(th) * v[:, None]))
    l = np.degrees(np.arctan2(pts[1], pts[0])) % 360.0
    b = np.degrees(np.arcsin(np.clip(pts[2], -1.0, 1.0)))
    return l, b


def main():
    results = []
    for name, fname, color in CATALOGUES:
        l, b, spin = gs.load_spin_catalogue(os.path.join(DATA, fname))
        R = gs.analyse(l, b, spin, cmb_axis_lb=CMB_AXIS)
        results.append((name, color, l, b, R))

    lines = ["Shared-axis test on two independent galaxy-spin catalogues",
             "=" * 60]
    for name, color, l, b, R in results:
        lines += [
            f"\n{name}",
            f"  CW fraction      = {R['f_cw']:.4f}  (excess {R['excess_percent']:+.2f}%,"
            f" binomial p = {R['binom_p']:.1e})",
            f"  dipole axis (l,b)= ({R['dipole_lb'][0]:.0f}, {R['dipole_lb'][1]:.0f})"
            f" deg  (bootstrap 68% = {R['axis_68_deg']:.0f} deg)",
            f"  dipole null p    = {R['dipole_null_p']:.3f}  (significance of a dipole"
            " BEYOND the monopole)",
            f"  angle to CMB axis= {R['angle_to_cmb_axis']:.0f} deg ;"
            f"  to Galactic pole = {R['angle_to_galactic_pole']:.0f} deg",
        ]
    # angle between the two datasets' axes
    from cartasis_sims.galaxy_spins import lb_to_vec, _acute_deg
    a0 = lb_to_vec(*results[0][4]["dipole_lb"])
    a1 = lb_to_vec(*results[1][4]["dipole_lb"])
    inter = _acute_deg(a0, a1)
    lines += [
        "",
        f"  GZ1 vs Iye axis separation = {inter:.0f} deg (both near the Galactic pole)",
        "",
        "READING: the two independent methods disagree on the asymmetry's sign and",
        "significance, neither yields a significant sky dipole, and both axes track the",
        "Galactic pole (~7 deg) far more than the CMB axis (~30 deg). Today's spin data",
        "is systematics-limited -- it neither confirms nor refutes the shared axis (a",
        "low-spin parent predicts a quiet sky too). The decisive test needs an all-sky,",
        "mirror-validated catalogue (Euclid, Rubin).",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "spin_datasets.txt"), "w") as f:
        f.write(text + "\n")

    # ---- figure ----------------------------------------------------------
    fig = plt.figure(figsize=(13.5, 5.6))
    axM = fig.add_subplot(1, 2, 1, projection="aitoff")
    axB = fig.add_subplot(1, 2, 2)

    axM.grid(True, alpha=0.3)
    # reference directions
    for (l, bb), lab, mk, c in [
            (CMB_AXIS, "CMB axis of evil", "*", "k"),
            ((0.0, 90.0), "N/S Galactic pole", "P", "0.5"),
            (ECL_POLE, "ecliptic pole", "D", "0.6")]:
        for sgn in (+1, -1):
            x, y = _aitoff_xy(l + (180 if sgn < 0 else 0), sgn * bb)
            axM.plot(x, y, mk, color=c, ms=12 if mk == "*" else 8,
                     label=lab if sgn > 0 else None)
    # the fitted dataset axes (headless: plot point and antipode) + bootstrap radius
    for name, color, l, b, R in results:
        la, ba = R["dipole_lb"]
        rad = R["axis_68_deg"]
        for sgn in (+1, -1):
            ll = la + (180 if sgn < 0 else 0)
            x, y = _aitoff_xy(ll, sgn * ba)
            axM.plot(x, y, "o", color=color, ms=9,
                     label=name.split(" (")[0] if sgn > 0 else None)
            # proper spherical bootstrap small-circle (broken at the l-wrap)
            cl, cb = _small_circle(ll, sgn * ba, rad)
            cx, cy = _aitoff_xy(cl, cb)
            cx[np.abs(np.diff(cx, prepend=cx[0])) > np.pi] = np.nan
            axM.plot(cx, cy, "-", color=color, lw=0.8, alpha=0.5)
    axM.set_title("Inferred spin axes vs reference directions\n"
                  "(Galactic, headless axes shown both ways)", fontsize=10.5)
    axM.legend(loc="upper right", bbox_to_anchor=(1.18, 1.18), fontsize=7.5)

    # right: asymmetry + dipole significance comparison
    ys = np.arange(len(results))[::-1]
    for y, (name, color, l, b, R) in zip(ys, results):
        n = R["n"]; f = R["f_cw"]
        err = 100.0 * np.sqrt(f * (1 - f) / n)
        axB.errorbar(R["excess_percent"], y, xerr=2 * err, fmt="o", color=color,
                     ms=8, capsize=4, lw=1.6)
        axB.text(R["excess_percent"], y + 0.16, name.split(" (")[0],
                 ha="center", fontsize=9, color=color)
        axB.text(R["excess_percent"], y - 0.22,
                 f"dipole null $p$={R['dipole_null_p']:.2f},  "
                 f"{R['angle_to_cmb_axis']:.0f}$^\\circ$ from CMB axis",
                 ha="center", fontsize=7.5, color="0.35")
    axB.axvline(0.0, color="k", ls="--", lw=1.2)
    axB.text(0.1, len(ys) - 0.5, "isotropy", fontsize=8.5)
    axB.set_xlim(-6, 6)
    axB.set_ylim(-0.6, len(ys) - 0.1)
    axB.set_yticks([])
    axB.set_xlabel("CW$-$CCW excess  [%]  ($\\pm2\\sigma$ binomial)")
    axB.set_title("Two methods, two answers\n(the systematic, not the cosmos)",
                  fontsize=10.5)
    axB.grid(True, axis="x", alpha=0.2)

    fig.suptitle("The shared-axis test on real, independent spin catalogues: "
                 "systematics-limited", fontsize=12.5, y=1.00)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(PDF_DIR, "spin_datasets.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "spin_datasets.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'spin_datasets.pdf')}")


if __name__ == "__main__":
    main()
