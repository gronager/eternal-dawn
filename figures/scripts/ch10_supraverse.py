#!/usr/bin/env python3
r"""The supraverse, gravity-scaled and polygonal: a still frame of the foam.

The mature foam is a Johnson-Mehl tessellation (void_foam.py): OGUs nucleate, grow
at the causal speed, and impinge into curved polygons. Rendered in grayscale -- no
colour, only void and eternal structure: each cell a gravity-scaled well (bright
core, dark rim), the Cartasis membranes between them toned by chirality (dark for
matter lineages, light for antimatter, inherited), with a hint of nested BHU
children in the larger cells. Our origin, frozen for the page.

Renders figures/pdf/supraverse.pdf and a high-res PNG.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt
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

    # per-cell properties: chirality (inherited per lineage), tones the membrane
    matter = rng.random(N) > 0.5

    # normalise the radial distance within each cell -> gravity-scaled well
    cell_max = np.ones(N)
    for i in range(N):
        m = labels == i
        if m.any():
            cell_max[i] = max(dmin[m].max(), 1e-3)
    r = dmin / cell_max[labels]               # 0 at core, 1 at rim

    # grayscale value: a luminous gravity well, bright core to dark rim,
    # with a faint per-cell tonal jitter so the foam breathes
    tone = rng.uniform(-0.05, 0.05, N)[labels]
    V = 0.82 - 0.54 * r ** 1.3 + tone
    V = np.clip(V, 0.12, 0.95)

    # Cartasis membranes: tone the cell boundaries by chirality
    # (dark seam for matter lineages, light seam for antimatter)
    edge = np.zeros(labels.shape, dtype=bool)
    edge[:-1, :] |= labels[:-1, :] != labels[1:, :]
    edge[1:, :] |= labels[:-1, :] != labels[1:, :]
    edge[:, :-1] |= labels[:, :-1] != labels[:, 1:]
    edge[:, 1:] |= labels[:, :-1] != labels[:, 1:]
    membrane_tone = np.where(matter[labels], 0.04, 0.96)
    V = np.where(edge, membrane_tone, V)

    fig = plt.figure(figsize=(16, 9), facecolor="black")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(V, origin="lower", extent=(0, aspect, 0, 1), cmap="gray",
              vmin=0.0, vmax=1.0, interpolation="bilinear")
    ax.set_xlim(0, aspect)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])

    # a hint of nesting: BHU children as faint rings in the larger cells
    order = np.argsort(-cell_max)
    for i in order[:14]:
        cx, cy, cr = xs[i], ys[i], cell_max[i]
        ring = "0.05" if matter[i] else "0.97"
        for _ in range(rng.integers(2, 5)):
            rho = cr * 0.55 * np.sqrt(rng.random())
            th = rng.random() * 2 * np.pi
            br = cr * rng.uniform(0.05, 0.12)
            ax.add_patch(Circle((cx + rho * np.cos(th), cy + rho * np.sin(th)), br,
                                facecolor="none", edgecolor=ring, lw=0.5, alpha=0.4))

    fig.text(0.022, 0.045, "THE SUPRAVERSE", color="white", fontsize=26,
             fontweight="bold", alpha=0.92, family="sans-serif")
    fig.text(0.022, 0.022,
             "gravity-scaled polygonal foam  ·  Johnson-Mehl cells (OGUs grown to "
             "impingement)  ·  dark seams = matter, light seams = antimatter "
             "lineages  ·  the seams are Cartasis  ·  no colour, only void and "
             "eternal structure",
             color="0.82", fontsize=9.5, alpha=0.9, family="sans-serif")

    out = os.path.join(PDF_DIR, "supraverse.pdf")
    fig.savefig(out, facecolor="black")
    fig.savefig(out.replace(".pdf", ".png"), dpi=170, facecolor="black")
    plt.close(fig)
    print(f"wrote {out} and .png  (cells={N}, grid={grid})")


if __name__ == "__main__":
    main()
