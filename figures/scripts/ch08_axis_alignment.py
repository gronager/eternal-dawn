#!/usr/bin/env python3
r"""Corner B: do the large-angle anomaly axes share one axis (a parent spin)?

Computes pairwise alignments and an orientation-tensor concentration test
(Monte-Carlo calibrated), renders a Galactic-coordinate sky map, and -- crucially
-- contrasts the cosmological reading (galaxy-spin axis tracks the CMB axis) with
the systematic reading (it tracks the Galactic pole).

Run with the sims venv, or via `make figures`.
"""
from __future__ import annotations

import json
import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import axes as ax

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

BY_NAME = {a.name: a for a in ax.AXES}
VEC = {a.name: ax.lb_to_vec(a.l, a.b) for a in ax.AXES}


def pair(n1: str, n2: str) -> dict:
    theta = ax.acute_angle_deg(VEC[n1], VEC[n2])
    return {"a": n1, "b": n2, "acute_deg": theta, "p_iso": ax.axis_pvalue(theta)}


def main() -> None:
    out: dict = {"pairwise": [], "notes": []}

    # --- headline comparisons -------------------------------------------------
    crux = [
        pair("Galaxy-spin asymmetry", "CMB axis of evil (l=2/3)"),
        pair("Galaxy-spin asymmetry", "Galactic pole"),
        pair("CMB axis of evil (l=2/3)", "CMB kinematic dipole"),
        pair("CMB axis of evil (l=2/3)", "North ecliptic pole"),
    ]
    out["pairwise"] = crux

    # --- concentration test on the well-sourced axes only ---------------------
    # (quadrupole/octupole are excluded: their per-multipole coordinates here are
    #  unverified placeholders, and the combined "axis of evil" already encodes
    #  their joint direction.)
    cosmo_names = ["CMB axis of evil (l=2/3)", "CMB kinematic dipole",
                   "CMB cold spot", "Galaxy-spin asymmetry"]
    cosmo_vecs = np.array([VEC[n] for n in cosmo_names])
    mc = ax.monte_carlo_pvalue(cosmo_vecs, n_trials=200_000, seed=1)
    pax = mc["principal_axis_vec"]
    pax_lb = (np.degrees(np.arctan2(pax[1], pax[0])) % 360,
              np.degrees(np.arcsin(np.clip(pax[2], -1, 1))))
    out["concentration"] = {
        "members": cosmo_names, "tau1": mc["tau1"], "p_value": mc["p_value"],
        "principal_axis_lb": [round(pax_lb[0], 1), round(abs(pax_lb[1]), 1)],
    }

    # --- print ----------------------------------------------------------------
    lines = ["Corner B: axis alignment of large-angle anomalies",
             "=" * 64,
             "Pairwise acute angles (isotropic p = 1 - cos theta):"]
    for p in crux:
        lines.append(f"  {p['a']:>26s}  vs  {p['b']:<26s}"
                     f"  {p['acute_deg']:5.1f} deg   p_iso={p['p_iso']:.3f}")
    lines += [
        "",
        f"Concentration of {len(cosmo_names)} CMB+galaxy axes:",
        f"  tau1 = {mc['tau1']:.3f}  (1/3 isotropic, 1 = perfectly aligned)",
        f"  Monte-Carlo p (isotropic) = {mc['p_value']:.4f}",
        f"  principal axis (l,b) = ({pax_lb[0]:.1f}, {abs(pax_lb[1]):.1f}) deg",
        "",
        "READING:",
        "  The galaxy-spin axis sits ~{:.0f} deg from the CMB axis but ~{:.0f} deg"
        .format(crux[0]['acute_deg'], crux[1]['acute_deg']),
        "  from the Galactic pole. Closer to the Galactic pole => the mundane",
        "  Milky-Way classification systematic is not excluded; the parent-spin",
        "  reading needs the galaxy axis to favour the CMB axis instead.",
        "  All p-values are a-posteriori-inflated lower bounds (selected axes).",
    ]
    table = "\n".join(lines)
    print(table)
    with open(os.path.join(OUT_DIR, "axis_alignment.txt"), "w") as f:
        f.write(table + "\n")
    with open(os.path.join(OUT_DIR, "axis_alignment.json"), "w") as f:
        json.dump(out, f, indent=2)

    # --- sky map --------------------------------------------------------------
    fig = plt.figure(figsize=(10, 5.6))
    axm = fig.add_subplot(111, projection="mollweide")
    colors = {"cmb": "C3", "galaxy": "C0", "geometry": "0.55"}
    markers = {"cmb": "o", "galaxy": "*", "geometry": "s"}
    seen = set()
    for a in ax.AXES:
        v = VEC[a.name]
        for sgn in (+1, -1):  # plot the axis and its antipode
            l = np.degrees(np.arctan2(sgn * v[1], sgn * v[0]))
            b = np.degrees(np.arcsin(np.clip(sgn * v[2], -1, 1)))
            lab = a.kind if a.kind not in seen else None
            seen.add(a.kind)
            axm.scatter(np.radians(((l + 180) % 360) - 180), np.radians(b),
                        c=colors[a.kind], marker=markers[a.kind],
                        s=90 if a.kind != "geometry" else 55,
                        edgecolors="k", linewidths=0.4, label=lab, zorder=3)
        if a.kind != "geometry":
            l0 = ((a.l + 180) % 360) - 180
            axm.annotate(a.name.replace("CMB ", "").replace(" asymmetry", ""),
                         (np.radians(l0), np.radians(a.b)),
                         textcoords="offset points", xytext=(5, 5), fontsize=7)
    # principal axis of the concentration test
    for sgn in (+1, -1):
        l = np.degrees(np.arctan2(sgn * pax[1], sgn * pax[0]))
        b = np.degrees(np.arcsin(np.clip(sgn * pax[2], -1, 1)))
        axm.scatter(np.radians(((l + 180) % 360) - 180), np.radians(b),
                    facecolors="none", edgecolors="C2", s=220, linewidths=1.8,
                    zorder=2, label="principal axis" if sgn == 1 else None)
    axm.grid(True, alpha=0.3)
    fig.suptitle("Large-angle anomaly axes in Galactic coordinates\n"
                 fr"concentration $\tau_1={mc['tau1']:.2f}$, "
                 fr"MC $p={mc['p_value']:.3f}$ (a-posteriori inflated)",
                 fontsize=11, y=0.99)
    axm.legend(loc="lower center", bbox_to_anchor=(0.5, -0.20), ncol=4,
               fontsize=8, frameon=False)
    fig.subplots_adjust(top=0.84, bottom=0.12)
    out_pdf = os.path.join(PDF_DIR, "axis_alignment.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
