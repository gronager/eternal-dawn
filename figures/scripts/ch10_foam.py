#!/usr/bin/env python3
r"""Foam wallpaper: the gravity-scaled condensed void, zooming to 'home'.

A gravity-scaled conformal rendering (Chapter 10): in the condensed void every
universe-bubble is drawn at comparable visual size regardless of physical scale.
We zoom: the void full of OG universes -> inside one OGU, a foam of black-hole
universes -> inside one BHU, sub-universes, one of which holds a pinprick: home.
Chirality is colour-coded (warm = matter, cool = antimatter).

Outputs figures/pdf/foam_wallpaper.pdf and a high-res PNG suitable as wallpaper.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, ConnectionPatch

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
os.makedirs(PDF_DIR, exist_ok=True)

BG = "#080810"
WARM = ["#ff7a33", "#ff9d3c", "#ffb347", "#e8602c", "#ffc66e"]   # matter
COOL = ["#3aa0ff", "#5bb8ff", "#7d6cff", "#46d6e0", "#8a7dff"]   # antimatter


def pack(n, rng, rmin, rmax, xlim=(0, 1), ylim=(0, 1), tries=600):
    """Random non-overlapping-ish circles (gravity-scaled: comparable radii)."""
    circ = []
    for _ in range(tries):
        if len(circ) >= n:
            break
        r = rng.uniform(rmin, rmax)
        x = rng.uniform(xlim[0] + r, xlim[1] - r)
        y = rng.uniform(ylim[0] + r, ylim[1] - r)
        if all((x - cx) ** 2 + (y - cy) ** 2 > (1.05 * (r + cr)) ** 2
               for cx, cy, cr, _ in circ):
            circ.append((x, y, r, rng.random() > 0.5))   # last: matter?
    return circ


def draw_foam(ax, circ, rng, vacuum_dots=400, edge="#cdd6ff"):
    ax.set_facecolor(BG)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    # faint quantum-vacuum fluctuation speckle
    ax.scatter(rng.random(vacuum_dots), rng.random(vacuum_dots),
               s=rng.uniform(0.2, 2.0, vacuum_dots), c="#9aa0c0",
               alpha=0.18, linewidths=0, zorder=0)
    for x, y, r, matter in circ:
        pal = WARM if matter else COOL
        col = pal[rng.integers(len(pal))]
        # soft halo (gravity glow) + membrane ring
        ax.add_patch(Circle((x, y), r, facecolor=col, alpha=0.16, lw=0, zorder=1))
        ax.add_patch(Circle((x, y), r, facecolor=col, alpha=0.30,
                            edgecolor=edge, lw=0.6, zorder=2))
    for spine in ax.spines.values():
        spine.set_color("#30364f")


def main() -> None:
    rng = np.random.default_rng(7)

    fig = plt.figure(figsize=(16, 9), facecolor=BG)
    axV = fig.add_axes([0.0, 0.0, 1.0, 1.0])          # the void (full bleed)
    void = pack(46, rng, 0.045, 0.085)
    draw_foam(axV, void, rng, vacuum_dots=900)

    # choose one OG universe to open
    ogu = void[len(void) // 2]
    ox, oy, orad = ogu[0], ogu[1], ogu[2]
    axV.add_patch(Circle((ox, oy), orad, facecolor="none",
                         edgecolor="white", lw=1.6, zorder=4))
    axV.text(ox, oy - orad - 0.018, "OG universe", color="white", ha="center",
             va="top", fontsize=9, alpha=0.9)

    # inset A: inside the OGU -> a foam of black-hole universes
    axA = fig.add_axes([0.60, 0.50, 0.30, 0.40])
    bhu_foam = pack(28, rng, 0.06, 0.11)
    draw_foam(axA, bhu_foam, rng, vacuum_dots=250, edge="#e7ecff")
    bhu = bhu_foam[len(bhu_foam) // 2 + 2]
    bx, by, brad = bhu[0], bhu[1], bhu[2]
    axA.add_patch(Circle((bx, by), brad, facecolor="none",
                         edgecolor="white", lw=1.4, zorder=4))
    axA.set_title("inside an OGU: black-hole universes", color="white",
                  fontsize=9, pad=4)

    # inset B: inside the BHU -> sub-universes, one is 'home'
    axB = fig.add_axes([0.60, 0.06, 0.30, 0.36])
    sub = pack(16, rng, 0.07, 0.13)
    draw_foam(axB, sub, rng, vacuum_dots=140, edge="#e7ecff")
    home = sub[len(sub) // 2]
    hx, hy, hrad, _ = home
    axB.add_patch(Circle((hx, hy), hrad, facecolor="none", edgecolor="white",
                         lw=1.4, zorder=4))
    # the pin
    axB.plot([hx], [hy], marker="o", ms=4, color="white", zorder=6)
    axB.annotate("home", (hx, hy), textcoords="offset points", xytext=(14, 10),
                 color="white", fontsize=11, fontweight="bold", zorder=6,
                 arrowprops=dict(arrowstyle="-", color="white", lw=1.0))
    axB.set_title("inside a BHU: a pinprick called home", color="white",
                  fontsize=9, pad=4)

    # callout connectors void->A->B
    for (xyA, axa, xyB, axb) in [((ox + orad * 0.7, oy + orad * 0.7), axV,
                                  (0.02, 0.98), axA),
                                 ((bx + brad * 0.7, by - brad * 0.7), axA,
                                  (0.02, 0.98), axB)]:
        con = ConnectionPatch(xyA=xyA, coordsA=axa.transData,
                              xyB=xyB, coordsB=axb.transAxes,
                              color="white", lw=0.7, alpha=0.5,
                              linestyle=(0, (4, 3)))
        fig.add_artist(con)

    fig.text(0.03, 0.06, "THE CONDENSED VOID", color="white", fontsize=22,
             fontweight="bold", alpha=0.92, family="sans-serif")
    fig.text(0.03, 0.03,
             "gravity-scaled foam  ·  warm = matter, cool = antimatter  "
             "·  home: inside a BHU, inside an OGU",
             color="#aab0d0", fontsize=11)

    out = os.path.join(PDF_DIR, "foam_wallpaper.pdf")
    fig.savefig(out, facecolor=BG)
    fig.savefig(out.replace(".pdf", ".png"), dpi=170, facecolor=BG)
    plt.close(fig)
    print(f"wrote {out} and .png")


if __name__ == "__main__":
    main()
