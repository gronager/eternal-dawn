#!/usr/bin/env python3
r"""The supraverse, gravity-scaled and polygonal: a still frame of the lava lamp.

The mature foam is a Johnson-Mehl tessellation (void_foam.py): OGUs nucleate, grow
at the causal speed, and impinge into curved polygons. Here it is rendered in colour
as the "condensed void" -- chirality split into two hue families (matter warm,
antimatter cool, inherited per lineage), each cell a gravity-scaled well (bright
core, dark rim) with glowing Cartasis membranes between, and a hint of nesting
(BHU children) in the larger cells. This is our origin, frozen for the page; the
same construction animated is the supraverse screensaver.

Renders figures/pdf/supraverse.pdf and a high-res PNG.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
from matplotlib.patches import Circle

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
os.makedirs(PDF_DIR, exist_ok=True)


def johnson_mehl(xs, ys, ts, grid, aspect=1.0):
    """Incremental Johnson-Mehl: returns (labels, dmin) on a grid x in [0,aspect],
    y in [0,1]. labels[p] = first seed to reach p; dmin[p] = its distance."""
    gx = (np.arange(int(grid * aspect)) + 0.5) / grid
    gy = (np.arange(grid) + 0.5) / grid
    X, Y = np.meshgrid(gx, gy)
    best = np.full(X.shape, np.inf)
    labels = np.zeros(X.shape, dtype=int)
    dmin = np.zeros(X.shape)
    for i in range(len(xs)):
        d = np.hypot(X - xs[i], Y - ys[i])
        a = ts[i] + d
        upd = a < best
        best = np.where(upd, a, best)
        labels = np.where(upd, i, labels)
        dmin = np.where(upd, d, dmin)
    return labels, dmin


def main() -> None:
    rng = np.random.default_rng(20)
    aspect = 16.0 / 9.0
    N = 120
    xs = rng.uniform(0, aspect, N)
    ys = rng.uniform(0, 1, N)
    ts = rng.uniform(0, 1, N) ** 1.4          # spread of birth times (older = larger)
    grid = 760

    labels, dmin = johnson_mehl(xs, ys, ts, grid, aspect)

    # per-cell properties
    matter = rng.random(N) > 0.5              # chirality, inherited per lineage
    hue_jit = rng.uniform(-0.035, 0.035, N)
    # matter: warm (red->gold ~0.02-0.11); antimatter: cool (teal->violet ~0.50-0.62)
    base_hue = np.where(matter, rng.uniform(0.01, 0.11, N),
                        rng.uniform(0.50, 0.64, N)) + hue_jit
    sat = np.where(matter, rng.uniform(0.72, 0.92, N), rng.uniform(0.55, 0.78, N))

    # normalise the radial distance within each cell -> gravity-scaled well
    cell_max = np.ones(N)
    for i in range(N):
        m = labels == i
        if m.any():
            cell_max[i] = max(dmin[m].max(), 1e-3)
    r = dmin / cell_max[labels]               # 0 at core, 1 at rim

    H = base_hue[labels]
    S = sat[labels] * (0.85 + 0.15 * r)       # a touch more saturated at the rim
    V = 0.96 - 0.62 * r ** 1.3                # bright core, dark rim (the well)

    # glowing Cartasis membranes: darken cell boundaries, then a thin bright lip
    edge = np.zeros(labels.shape, dtype=bool)
    edge[:-1, :] |= labels[:-1, :] != labels[1:, :]
    edge[1:, :] |= labels[:-1, :] != labels[1:, :]
    edge[:, :-1] |= labels[:, :-1] != labels[:, 1:]
    edge[:, 1:] |= labels[:, :-1] != labels[:, 1:]
    V = np.where(edge, 0.06, V)               # dark membrane gap

    rgb = hsv_to_rgb(np.dstack([H % 1.0, np.clip(S, 0, 1), np.clip(V, 0, 1)]))

    fig = plt.figure(figsize=(16, 9), facecolor="black")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(rgb, origin="lower", extent=(0, aspect, 0, 1), interpolation="bilinear")
    ax.set_xlim(0, aspect)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])

    # a hint of nesting: BHU children as faint glowing rings in the larger cells
    order = np.argsort(-cell_max)
    for i in order[:14]:
        cx, cy, cr = xs[i], ys[i], cell_max[i]
        ring = "1.0" if matter[i] else "0.85"
        for k in range(rng.integers(2, 5)):
            rho = cr * 0.55 * np.sqrt(rng.random())
            th = rng.random() * 2 * np.pi
            br = cr * rng.uniform(0.05, 0.12)
            ax.add_patch(Circle((cx + rho * np.cos(th), cy + rho * np.sin(th)), br,
                                facecolor="none", edgecolor=ring, lw=0.5, alpha=0.35))

    fig.text(0.022, 0.045, "THE SUPRAVERSE", color="white", fontsize=26,
             fontweight="bold", alpha=0.92, family="sans-serif")
    fig.text(0.022, 0.022,
             "gravity-scaled polygonal foam  ·  Johnson-Mehl cells (OGUs grown to "
             "impingement)  ·  warm = matter, cool = antimatter lineages  ·  glowing "
             "membranes are Cartasis  ·  our origin, one lit pixel among the "
             "uncountably many",
             color="0.85", fontsize=9.5, alpha=0.9, family="sans-serif")

    out = os.path.join(PDF_DIR, "supraverse.pdf")
    fig.savefig(out, facecolor="black")
    fig.savefig(out.replace(".pdf", ".png"), dpi=170, facecolor="black")
    plt.close(fig)
    print(f"wrote {out} and .png  (cells={N}, grid={grid})")


if __name__ == "__main__":
    main()
