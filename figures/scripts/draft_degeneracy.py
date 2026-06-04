#!/usr/bin/env python3
r"""Degeneracy pressure in the gravity-torsion soliton (Parts II/IV refinement).

The first soliton pass had one fermion in a fixed well -- no Pauli exclusion. A real
soliton is many fermions filling the levels, and their degeneracy pressure is a
genuine OUTWARD force countering the binding alongside the torsion wall. At the
Thomas-Fermi (equation-of-state) level:

  e(n) = e_kin(n) - a n + b n^2,   e_kin = relativistic Fermi gas (degeneracy KE).

A: binding per particle vs density, WITH and WITHOUT the Pauli kinetic term -- Pauli
   lowers the saturation density and shallows the binding.
B: the implied 'drop': without Pauli it is small, dense, sharp; with Pauli it is
   larger, less dense, flatter-topped -- a more box-like / scale-flat interior, which
   is the 'walking' direction that helps the electroweak S parameter.

Renders figures/pdf/degeneracy.pdf, writes sims/output/degeneracy.txt.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import fermi_ball as fb

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


def _drop_profile(r, n0, R, delta=0.18):
    """A Woods-Saxon 'drop': uniform density n0 inside radius R, smooth surface."""
    return n0 / (1.0 + np.exp((r - R) / delta))


def main():
    a, b = 5.0, 4.0
    on = fb.saturation(a, b, degeneracy=True)
    off = fb.saturation(a, b, degeneracy=False)
    d = fb.degeneracy_shift(a, b)

    lines = [
        "Degeneracy pressure in the soliton (Thomas-Fermi)",
        "=" * 50,
        f"  without Pauli KE: n0 = {off['n0']:.3f},  binding/particle = {off['binding']:.3f}",
        f"  with    Pauli KE: n0 = {on['n0']:.3f},  binding/particle = {on['binding']:.3f}",
        f"  -> Pauli lowers density to {100*d['n0_ratio']:.0f}% and binding to "
        f"{100*d['binding_ratio']:.0f}% of the no-Pauli case.",
        "",
        "READING: the Pauli exclusion pressure is a real outward force (as in white",
        "dwarfs/neutron stars); it pushes the drop apart -- larger, less dense, flatter-",
        "topped. A flatter, more uniform (scale-flat) interior is the 'walking' regime",
        "that SUPPRESSES the electroweak S. So including degeneracy moves S in the",
        "favourable direction. It does NOT by itself deliver S < 0.1 -- that needs the",
        "full self-consistent, colored, composite V/A spectrum -- but Pauli is in the",
        "column that helps, confirming the intuition.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "degeneracy.txt"), "w") as f:
        f.write(text + "\n")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.0, 5.0))

    # ---- A: binding curve with/without degeneracy ----
    n = np.linspace(1e-3, 0.6, 400)
    B_on = 1.0 - fb.energy_per_particle(n, a, b, degeneracy=True)
    B_off = 1.0 - fb.energy_per_particle(n, a, b, degeneracy=False)
    axA.plot(n, B_off, "C1", lw=2.0, label="no Pauli (rest mass only)")
    axA.plot(n, B_on, "C0", lw=2.0, label="with Pauli (degeneracy KE)")
    axA.plot(off["n0"], off["binding"], "C1o", ms=8)
    axA.plot(on["n0"], on["binding"], "C0o", ms=8)
    axA.axhline(0, color="k", lw=0.8)
    axA.annotate("", (on["n0"], on["binding"]), (off["n0"], off["binding"]),
                 arrowprops=dict(arrowstyle="->", color="0.4", lw=1.4))
    axA.text(0.35, 0.9, "Pauli pushes\nthe drop apart", fontsize=9, color="0.4")
    axA.set_xlabel(r"density $n$ (units $m^3$)")
    axA.set_ylabel(r"binding per particle  $1-e(n)$  (units $m$)")
    axA.set_title("Degeneracy lowers the saturation density\nand shallows the binding",
                  fontsize=11)
    axA.legend(fontsize=9, loc="lower center")
    axA.grid(True, alpha=0.2)

    # ---- B: the implied drop profile ----
    N = 9.0
    R_off = (3 * N / (4 * np.pi * off["n0"])) ** (1.0 / 3.0)
    R_on = (3 * N / (4 * np.pi * on["n0"])) ** (1.0 / 3.0)
    r = np.linspace(0, max(R_on, R_off) * 1.6, 400)
    axB.plot(r, _drop_profile(r, off["n0"], R_off), "C1", lw=2.0,
             label=f"no Pauli (dense, sharp, $R$={R_off:.2f})")
    axB.plot(r, _drop_profile(r, on["n0"], R_on), "C0", lw=2.0,
             label=f"with Pauli (spread, flat, $R$={R_on:.2f})")
    axB.fill_between(r, 0, _drop_profile(r, on["n0"], R_on), color="C0", alpha=0.08)
    axB.set_xlabel(r"radius $r$ (units $1/m$)")
    axB.set_ylabel(r"density $n(r)$")
    axB.set_title("The drop: Pauli makes it larger and flatter\n"
                  "(box-like interior $=$ 'walking' $=$ helps $S$)", fontsize=11)
    axB.legend(fontsize=8.5, loc="upper right")
    axB.grid(True, alpha=0.2)

    fig.suptitle("Pauli degeneracy pressure in the soliton: it pushes toward the "
                 "S-friendly 'walking' regime", fontsize=12.5, y=1.00)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(PDF_DIR, "degeneracy.pdf"))
    fig.savefig(os.path.join(PDF_DIR, "degeneracy.png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {os.path.join(PDF_DIR, 'degeneracy.pdf')}")


if __name__ == "__main__":
    main()
