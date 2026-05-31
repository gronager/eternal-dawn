#!/usr/bin/env python3
r"""De-confound the galaxy-spin axis on real Galaxy Zoo 1 data.

Fits the spin dipole ourselves, bootstraps its axis, tests it against a
label-shuffle null, compares it to the CMB axis and the Galactic pole, and
exposes the Galactic-latitude systematic. Renders figures/pdf/galaxy_dipole.pdf.
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
CSV = os.path.join(ROOT, "data", "galaxy_spins", "gz1_clean.csv")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def load_gz1(path):
    d = np.genfromtxt(path, delimiter=",", names=True)
    ra, dec, pcS, paS = d["RAJ2000"], d["DEJ2000"], d["pcS"], d["paS"]
    spin = np.where(pcS > paS, 1.0, -1.0)
    l, b = gs.equatorial_to_galactic(ra, dec)
    return np.asarray(l), np.asarray(b), spin


def main() -> None:
    l, b, spin = load_gz1(CSV)
    r = gs.analyse(l, b, spin)

    lines = [
        "De-confounding the galaxy-spin axis (Galaxy Zoo 1, clean spirals)",
        "=" * 66,
        f"  galaxies (CW / CCW)        = {r['n']}  ({r['n_cw']} / {r['n_ccw']})",
        f"  clockwise fraction         = {r['f_cw']:.4f}  "
        f"(excess {r['excess_percent']:+.2f}%)",
        f"  binomial p (asym != 0.5)   = {r['binom_p']:.2e}",
        "",
        f"  fitted dipole axis (l,b)   = ({r['dipole_lb'][0]:.1f}, "
        f"{r['dipole_lb'][1]:.1f}) deg",
        f"  bootstrap 68% axis radius  = {r['axis_68_deg']:.1f} deg",
        f"  dipole amplitude (3|D|)    = {r['amplitude_3D']:.4f}",
        f"  label-shuffle null p       = {r['dipole_null_p']:.3f}  "
        "(p<0.05 => a real dipole beyond the monopole)",
        "",
        f"  angle: dipole vs CMB axis      = {r['angle_to_cmb_axis']:.1f} deg",
        f"  angle: dipole vs Galactic pole = {r['angle_to_galactic_pole']:.1f} deg",
        "",
        "READING: Galaxy Zoo 1 is an SDSS (northern-cap) footprint with a known",
        "human handedness/perception bias (Land et al. 2008). Any asymmetry here",
        "is dominated by systematics, not cosmology; the dipole axis is mask- and",
        "bias-limited. A genuine test of the parent-spin prediction needs an",
        "all-sky, bias-controlled catalogue. This is the method, run on the",
        "cautionary real dataset.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "galaxy_dipole.txt"), "w") as f:
        f.write(text + "\n")

    fig = plt.figure(figsize=(11, 4.8))

    # Panel A: sky map (Galactic), galaxies colored by spin + key axes.
    axm = fig.add_subplot(121, projection="mollweide")
    rng = np.random.default_rng(0)
    sub = rng.choice(len(l), size=min(5000, len(l)), replace=False)
    ll = ((l[sub] + 180.0) % 360.0) - 180.0
    axm.scatter(np.radians(ll), np.radians(b[sub]),
                c=np.where(spin[sub] > 0, "C3", "C0"), s=2, alpha=0.35,
                linewidths=0)
    for lb, mk, col, lab in [(r["dipole_lb"], "*", "k", "spin dipole"),
                             ((260.0, 60.0), "o", "C2", "CMB axis"),
                             ((0.0, 90.0), "s", "C1", "Gal. pole")]:
        for sgn in (1, -1):
            v = gs.lb_to_vec(*lb) * sgn
            ll2 = np.degrees(np.arctan2(v[1], v[0]))
            bb2 = np.degrees(np.arcsin(np.clip(v[2], -1, 1)))
            axm.scatter(np.radians(((ll2 + 180) % 360) - 180), np.radians(bb2),
                        marker=mk, c=col, s=130, edgecolors="w", linewidths=0.6,
                        label=lab if sgn == 1 else None, zorder=5)
    axm.grid(True, alpha=0.3)
    axm.set_title("GZ1 clean spirals (red=CW, blue=CCW)\nGalactic coordinates",
                  fontsize=10)
    axm.legend(loc="lower center", bbox_to_anchor=(0.5, -0.22), ncol=3,
               fontsize=8, frameon=False)

    # Panel B: the systematic discriminant -- CW fraction vs |Galactic latitude|.
    ax2 = fig.add_subplot(122)
    absb = np.abs(b)
    edges = np.linspace(absb.min(), absb.max(), 9)
    cen, frac, err = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (absb >= lo) & (absb < hi)
        if m.sum() < 50:
            continue
        k = (spin[m] > 0).sum()
        nn = m.sum()
        cen.append(0.5 * (lo + hi))
        frac.append(k / nn)
        err.append(np.sqrt(0.25 / nn))
    ax2.errorbar(cen, frac, yerr=err, fmt="o-", color="C4", capsize=3)
    ax2.axhline(0.5, color="0.5", ls=":", lw=1.0)
    ax2.set_xlabel(r"$|b|$  Galactic latitude [deg]")
    ax2.set_ylabel("clockwise fraction")
    ax2.set_title("Discriminant: does the asymmetry track\nGalactic latitude?",
                  fontsize=10)
    ax2.grid(True, alpha=0.25)

    fig.suptitle("De-confounding the galaxy-spin axis (Galaxy Zoo 1, real data)",
                 fontsize=12, y=0.99)
    fig.subplots_adjust(left=0.04, right=0.96, wspace=0.22, top=0.84, bottom=0.14)
    out_pdf = os.path.join(PDF_DIR, "galaxy_dipole.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
