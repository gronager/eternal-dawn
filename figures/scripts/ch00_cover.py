#!/usr/bin/env python3
r"""Front cover: the recursive supraverse foam, portrait, full-bleed.

A clean #808080 gray supraverse, dusted with original-generation universes (OGUs)
whose interiors are the same gray and whose edges are black (matter) or white
(antimatter) in roughly equal numbers. Inside each OGU sit smaller BHU1s --- gray
inside, edged in the SAME colour as their parent --- and inside each BHU1 a few
dots: the BHU2s, the black holes we watch form in our own sky. Title in white.
Renders figures/pdf/cover.pdf (full-bleed) for the book cover page and the EPUB
cover image.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
os.makedirs(PDF_DIR, exist_ok=True)

GRAY = "#808080"       # the supraverse background AND every interior
BLACK = "0.04"         # matter-lineage edge
WHITE = "0.97"         # antimatter-lineage edge
TITLE_COLOR = "white"  # title text


def _pack(rng, R, lo, hi, n, gap, tries=400):
    """Pack up to n non-overlapping discs (radii in [lo, hi]) within radius R of
    the origin, by rejection. Returns a list of (dx, dy, rad)."""
    out = []
    for _ in range(n):
        for _ in range(tries):
            rad = rng.uniform(lo, hi)
            rho = (R - rad) * np.sqrt(rng.random())
            th = rng.random() * 2 * np.pi
            dx, dy = rho * np.cos(th), rho * np.sin(th)
            if all((dx - qx) ** 2 + (dy - qy) ** 2 > (rad + qr + gap) ** 2
                   for qx, qy, qr in out):
                out.append((dx, dy, rad))
                break
    return out


def draw_ogu(ax, x, y, r, color, rng):
    """An OGU: a gray disc edged in `color` (black=matter / white=antimatter),
    holding NON-overlapping same-coloured BHU1s, each holding NON-overlapping
    BHU2 dots."""
    ax.add_patch(Ellipse((x, y), 2 * r, 2 * r, facecolor=GRAY,
                         edgecolor=color, lw=1.7, zorder=3))
    # BHU1s inside, packed without overlap, same lineage colour as the parent
    for dx, dy, br in _pack(rng, 0.86 * r, 0.16 * r, 0.24 * r,
                            n=int(rng.integers(3, 6)), gap=0.03 * r):
        bx, by = x + dx, y + dy
        ax.add_patch(Ellipse((bx, by), 2 * br, 2 * br, facecolor=GRAY,
                             edgecolor=color, lw=0.9, zorder=4))
        # BHU2 dots inside each BHU1 (the black holes we observe), also non-overlapping
        for ex, ey, dr in _pack(rng, 0.80 * br, 0.13 * br, 0.20 * br,
                                n=int(rng.integers(2, 5)), gap=0.02 * br):
            ax.add_patch(Ellipse((bx + ex, by + ey), 2 * dr, 2 * dr,
                                facecolor=color, edgecolor="none", zorder=5))


def main() -> None:
    rng = np.random.default_rng(7)
    W, Hh = 8.27, 11.69                 # A4 portrait
    aspect = Hh / W
    fig = plt.figure(figsize=(W, Hh), facecolor=GRAY)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, aspect)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_facecolor(GRAY)
    for s in ax.spines.values():
        s.set_visible(False)

    # balanced black/white lineage colours, shuffled
    n_ogu = 12
    colors = [BLACK, WHITE] * (n_ogu // 2)
    rng.shuffle(colors)

    # sparse, well-separated OGUs; keep the central title band clear
    placed, ci = [], 0
    for _ in range(4000):
        if ci >= n_ogu:
            break
        r = rng.uniform(0.075, 0.115)
        x, y = rng.uniform(r, 1 - r), rng.uniform(r, aspect - r)
        if 0.50 * aspect < y < 0.80 * aspect and 0.12 < x < 0.88:
            continue                                  # leave room for the title
        if y < 0.13 * aspect and 0.12 < x < 0.88:
            continue                                  # leave room for the author line
        if all((x - px) ** 2 + (y - py) ** 2 > (1.9 * (r + pr)) ** 2
               for px, py, pr in placed):
            placed.append((x, y, r))
            draw_ogu(ax, x, y, r, colors[ci],
                     np.random.default_rng(int(abs(x * 1e6 + y * 1e3))))
            ci += 1

    # --- title block, white ---
    ax.text(0.5, aspect * 0.71, "ETERNAL\nDAWN", color=TITLE_COLOR,
            ha="center", va="center", fontsize=58, fontweight="bold",
            family="serif", linespacing=1.05, zorder=10)
    ax.text(0.5, aspect * 0.565, "A continuous, conservative cosmology",
            color=TITLE_COLOR, ha="center", va="center", fontsize=18, style="italic",
            family="serif", alpha=0.95, zorder=10)
    ax.text(0.5, aspect * 0.085, "MICHAEL GRONAGER, PhD", color=TITLE_COLOR,
            ha="center", va="center", fontsize=16, fontweight="bold",
            family="serif", zorder=10)

    out = os.path.join(PDF_DIR, "cover.pdf")
    fig.savefig(out, facecolor=GRAY, dpi=300)
    fig.savefig(out.replace(".pdf", ".png"), dpi=200, facecolor=GRAY)
    plt.close(fig)
    print(f"wrote {out} and .png  (OGUs={len(placed)})")


if __name__ == "__main__":
    main()
