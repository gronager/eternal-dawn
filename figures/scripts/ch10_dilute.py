#!/usr/bin/env python3
r"""The dilute supraverse, true to scale: round islands in a weighted void.

The instanton verdict (I >> I_crit) puts us in the dilute regime: OGUs are isolated,
round, and astronomically far apart. The vacuum has weight (rho_Lambda ~ 6e-27 kg/m^3),
so the void occupies real gravity-scaled volume -- universes do not touch (not
polygons) -- and because the regime is dilute, that void overwhelmingly dominates.
Grayscale, muted: a vast textured void with sparse round island universes (gravity
wells, chirality-toned Cartasis membranes), and an inset zoom (the gravity-scaling
needed to see one at all). Renders figures/pdf/dilute_supraverse.pdf and a PNG.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, ConnectionPatch

from cartasis_sims import void_scale as vs

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
os.makedirs(PDF_DIR, exist_ok=True)

VOID = "0.30"        # the weighted void (mid-dark gray)
BLACK = "0.04"
WHITE = "0.97"


def draw_island(ax, x, y, r, matter, rng, rings=True):
    """A round island universe: a gravity well with a chirality-toned membrane."""
    n = 60
    for i in range(n, 0, -1):
        rr = r * i / n
        val = 0.30 + 0.50 * (1 - (i / n) ** 1.3)      # bright core, dark rim
        ax.add_patch(Circle((x, y), rr, facecolor=str(min(val, 0.9)),
                            edgecolor="none", zorder=2))
    edge = BLACK if matter else WHITE
    ax.add_patch(Circle((x, y), r, facecolor="none", edgecolor=edge, lw=1.1,
                        zorder=3))
    if rings:
        for _ in range(rng.integers(2, 4)):            # faint nested BHU children
            rho = r * 0.5 * np.sqrt(rng.random())
            th = rng.random() * 2 * np.pi
            ax.add_patch(Circle((x + rho * np.cos(th), y + rho * np.sin(th)),
                                r * rng.uniform(0.06, 0.13), facecolor="none",
                                edgecolor=edge, lw=0.5, alpha=0.4, zorder=4))


def main() -> None:
    rng = np.random.default_rng(7)
    aspect = 16.0 / 9.0

    fig = plt.figure(figsize=(16, 9), facecolor=VOID)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, aspect)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor(VOID)

    # the weighted void: faint zero-point-energy speckle (it has mass, hence texture)
    sp = rng.uniform(0, aspect, 2600), rng.uniform(0, 1, 2600)
    ax.scatter(*sp, s=rng.uniform(0.2, 1.4, 2600),
               c=[str(s) for s in rng.uniform(0.22, 0.40, 2600)], alpha=0.5,
               linewidths=0, zorder=1)

    # sparse round island universes (dilute -- far apart, small in the frame)
    n_isl = 9
    placed = []
    for _ in range(n_isl * 40):
        if len(placed) >= n_isl:
            break
        r = rng.uniform(0.018, 0.05)
        x, y = rng.uniform(r, aspect - r), rng.uniform(r, 1 - r)
        if all((x - px) ** 2 + (y - py) ** 2 > (6.0 * (r + pr)) ** 2
               for px, py, pr in placed):           # enforce wide separation
            placed.append((x, y, r))
            draw_island(ax, x, y, r, rng.random() > 0.5, rng)

    # inset: zoom into one island (gravity-scaling -- the only way to see structure)
    hx, hy, hr = max(placed, key=lambda p: p[2])
    axm = fig.add_axes([0.66, 0.05, 0.31, 0.40])
    axm.set_facecolor("0.27")
    axm.set_xlim(hx - hr * 1.8, hx + hr * 1.8)
    axm.set_ylim(hy - hr * 1.8, hy + hr * 1.8)
    axm.set_xticks([]); axm.set_yticks([])
    for s in axm.spines.values():
        s.set_color(BLACK); s.set_linewidth(1.0)
    draw_island(axm, hx, hy, hr, True, np.random.default_rng(1))
    axm.plot([hx], [hy], "o", ms=4, color=BLACK, markeredgecolor=WHITE,
             markeredgewidth=0.7, zorder=20)
    axm.annotate("home", (hx, hy), textcoords="offset points", xytext=(12, 12),
                 color=BLACK, fontsize=11, fontweight="bold", zorder=20,
                 arrowprops=dict(arrowstyle="-", color=BLACK, lw=1))
    axm.set_title("gravity-scaled zoom: one round island (a BHU1 in its OGU)",
                  color=BLACK, fontsize=9, pad=3)
    con = ConnectionPatch((hx, hy), (hx, hy), "data", "data", axesA=ax, axesB=axm,
                          color=BLACK, lw=0.8, alpha=0.5, ls=(0, (4, 3)))
    fig.add_artist(con)

    fig.text(0.022, 0.05, "THE DILUTE SUPRAVERSE", color=WHITE, fontsize=24,
             fontweight="bold", alpha=0.95)
    fig.text(0.022, 0.025,
             f"true to scale: the void weighs $\\rho_\\Lambda\\approx${vs.vacuum_density():.0e}"
             " kg/m$^3$ (~3 protons/m$^3$), so it occupies real space -- universes are "
             "ROUND and isolated, not polygons; the void dominates, and they sit "
             "$\\sim e^{I/4}$ horizons apart  ($I\\sim10^{84}$).",
             color="0.88", fontsize=9.5, alpha=0.9)

    out = os.path.join(PDF_DIR, "dilute_supraverse.pdf")
    fig.savefig(out, facecolor=VOID)
    fig.savefig(out.replace(".pdf", ".png"), dpi=170, facecolor=VOID)
    plt.close(fig)
    print(f"wrote {out} and .png  (islands={len(placed)})")


if __name__ == "__main__":
    main()
