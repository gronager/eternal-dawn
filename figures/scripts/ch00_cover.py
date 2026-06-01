#!/usr/bin/env python3
r"""Front cover: the dilute supraverse, portrait, full-bleed, Bordeaux text.

A portrait (A4-proportion) version of the gravity-scaled dilute foam (ch10_dilute),
sized to fill the entire front cover, with the title set in Bordeaux. Renders
figures/pdf/cover.pdf (full-bleed, no margins) for use as the book's cover page and as
the EPUB cover image.
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

VOID = "0.30"          # weighted-void gray
BLACK = "0.04"
WHITE = "0.97"
BORDEAUX = "#6e1023"   # deep Bordeaux for the title text


def draw_island(ax, x, y, r, matter, rng):
    """A flat, slightly-oblate island: gravity-scaled disc + chirality-toned membrane."""
    ecc = rng.uniform(0.0, 0.26)
    ang = rng.uniform(0.0, 180.0)
    w = 2 * r * (1.0 + 0.5 * ecc)
    h = 2 * r * (1.0 - 0.5 * ecc)
    fill = str(0.40 + rng.uniform(-0.04, 0.04))
    ax.add_patch(Ellipse((x, y), w, h, angle=ang, facecolor=fill, edgecolor="none",
                         zorder=2))
    edge = BLACK if matter else WHITE
    ax.add_patch(Ellipse((x, y), w, h, angle=ang, facecolor="none", edgecolor=edge,
                         lw=1.1, zorder=3))
    th = np.deg2rad(ang + 90.0)
    ax.plot([x - 0.5 * h * np.cos(th), x + 0.5 * h * np.cos(th)],
            [y - 0.5 * h * np.sin(th), y + 0.5 * h * np.sin(th)],
            color=edge, lw=0.4, alpha=0.4, zorder=3)
    for _ in range(rng.integers(2, 4)):
        rho = r * 0.45 * np.sqrt(rng.random())
        tth = rng.random() * 2 * np.pi
        ax.add_patch(Ellipse((x + rho * np.cos(tth), y + rho * np.sin(tth)),
                            w * rng.uniform(0.06, 0.12), h * rng.uniform(0.06, 0.12),
                            angle=ang, facecolor="none", edgecolor=edge, lw=0.4,
                            alpha=0.4, zorder=4))


def main() -> None:
    rng = np.random.default_rng(20)
    # A4 portrait proportion (210 x 297 mm) -> aspect h/w = 297/210
    W, Hh = 8.27, 11.69
    aspect = Hh / W                                    # y-range per unit x
    fig = plt.figure(figsize=(W, Hh), facecolor=VOID)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, aspect)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_facecolor(VOID)
    for s in ax.spines.values():
        s.set_visible(False)

    # weighted-void speckle (zero-point energy)
    sp = rng.uniform(0, 1, 2600), rng.uniform(0, aspect, 2600)
    ax.scatter(*sp, s=rng.uniform(0.2, 1.4, 2600),
               c=[str(s) for s in rng.uniform(0.22, 0.40, 2600)], alpha=0.5,
               linewidths=0, zorder=1)

    # sparse, well-separated islands (dilute; illustrative sizes)
    placed = []
    for _ in range(400):
        if len(placed) >= 11:
            break
        r = rng.uniform(0.022, 0.060)
        x, y = rng.uniform(r, 1 - r), rng.uniform(r, aspect - r)
        if all((x - px) ** 2 + (y - py) ** 2 > (5.0 * (r + pr)) ** 2
               for px, py, pr in placed):
            placed.append((x, y, r))
            draw_island(ax, x, y, r, rng.random() > 0.5, np.random.default_rng(
                int(abs(x * 1e6 + y * 1e3))))

    # --- title block, Bordeaux ---
    ax.text(0.5, aspect * 0.70, "DAWN OF\nETERNITY", color=BORDEAUX,
            ha="center", va="center", fontsize=58, fontweight="bold",
            family="serif", linespacing=1.05, zorder=10)
    ax.text(0.5, aspect * 0.55, "Supraverse Cartasis Theory", color=BORDEAUX,
            ha="center", va="center", fontsize=20, style="italic",
            family="serif", alpha=0.95, zorder=10)
    ax.text(0.5, aspect * 0.085, "MICHAEL GRONAGER", color=BORDEAUX,
            ha="center", va="center", fontsize=16, fontweight="bold",
            family="serif", zorder=10)

    out = os.path.join(PDF_DIR, "cover.pdf")
    fig.savefig(out, facecolor=VOID, dpi=300)
    fig.savefig(out.replace(".pdf", ".png"), dpi=200, facecolor=VOID)
    plt.close(fig)
    print(f"wrote {out} and .png  (islands={len(placed)})")


if __name__ == "__main__":
    main()
