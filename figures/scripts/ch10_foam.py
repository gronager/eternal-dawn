#!/usr/bin/env python3
r"""Foam wallpaper: the gravity-scaled condensed void with its statistics baked in.

A self-consistent illustration of the thesis (Chapters 5, 10):

* OG universes (OGUs) fill the void with a *defined-median* size distribution --
  births outrun runaway growth, so many small + few large (log-normal radii).
* Each OGU spins (a median viable spin sets the median asymmetry, hence C); the
  larger ones carry a faint spin-axis tick.
* Nesting is self-similar at diameter ratio ~1/13 (a child universe is 1/12-1/15
  the size of its parent on the log-gravity scale -- NOT the ~1e-12 mass ratio),
  only a couple of layers deep before sub-bubbles are too small to host holes.
* Every disc is a universe; its periphery is the Cartasis membrane. Colour is set
  per bounce (warm = matter, cool = antimatter), random sign -- not inherited,
  because the bounce is CPT-even.
* A pin marks 'home', inside a black-hole universe inside an OGU; a magnifier
  makes it legible.

Outputs figures/pdf/foam_wallpaper.pdf and a high-res PNG.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, ConnectionPatch

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
os.makedirs(PDF_DIR, exist_ok=True)

BG = "#070710"
WARM = ["#ff7a33", "#ff9d3c", "#ffb347", "#e8602c"]   # matter
COOL = ["#3aa0ff", "#5bb8ff", "#7d6cff", "#46d6e0"]   # antimatter
NEST = 1.0 / 13.0          # child / parent diameter ratio (log-gravity scale)
W, H = 1.78, 1.0           # 16:9 void canvas


@dataclass
class Node:
    x: float
    y: float
    r: float
    level: int
    matter: bool
    children: list = field(default_factory=list)


def _color(rng, matter):
    pal = WARM if matter else COOL
    return pal[rng.integers(len(pal))]


def grow_children(node, rng, max_level, n_for_level):
    """Place children inside `node` at radius NEST*node.r, no heavy overlap."""
    if node.level >= max_level:
        return
    child_r = node.r * NEST
    if child_r < 0.0006:
        return
    n = n_for_level(node.level, rng)
    placed = []
    for _ in range(n * 12):
        if len(placed) >= n:
            break
        rho = (node.r - child_r * 1.3) * np.sqrt(rng.random())
        th = rng.random() * 2 * np.pi
        cx, cy = node.x + rho * np.cos(th), node.y + rho * np.sin(th)
        if all((cx - p.x) ** 2 + (cy - p.y) ** 2 > (1.15 * (child_r + p.r)) ** 2
               for p in placed):
            ch = Node(cx, cy, child_r, node.level + 1, rng.random() > 0.5)
            placed.append(ch)
            grow_children(ch, rng, max_level, n_for_level)
    node.children = placed


def sample_ogus(rng, n=34, median_r=0.052, sigma=0.5):
    """Log-normal OGU radii (defined median), rejection-packed into the void."""
    radii = np.sort(rng.lognormal(mean=np.log(median_r), sigma=sigma, size=n*3))[::-1]
    nodes = []
    for r in radii:
        if len(nodes) >= n:
            break
        r = float(min(r, 0.16))
        for _ in range(140):
            x, y = rng.uniform(r, W - r), rng.uniform(r, H - r)
            if all((x - p.x) ** 2 + (y - p.y) ** 2 > (1.04 * (r + p.r)) ** 2
                   for p in nodes):
                nodes.append(Node(x, y, r, 0, rng.random() > 0.5))
                break
    return nodes


def draw_node(ax, node, rng, glow=True):
    col = _color(rng, node.matter)
    if glow and node.r > 0.02:
        ax.add_patch(Circle((node.x, node.y), node.r * 1.18, facecolor=col,
                            alpha=0.10, lw=0, zorder=node.level))
    ax.add_patch(Circle((node.x, node.y), node.r, facecolor=col,
                        alpha=0.22 if node.level == 0 else 0.34,
                        edgecolor="#e8ecff",
                        lw=0.8 if node.level == 0 else 0.5, zorder=node.level + 1))
    # spin-axis tick on larger OGUs (the median viable spin that sets C)
    if node.level == 0 and node.r > 0.045:
        ang = rng.random() * np.pi
        dx, dy = node.r * 0.7 * np.cos(ang), node.r * 0.7 * np.sin(ang)
        ax.plot([node.x - dx, node.x + dx], [node.y - dy, node.y + dy],
                color="#e8ecff", lw=0.6, alpha=0.5, zorder=node.level + 2)
    for ch in node.children:
        draw_node(ax, ch, rng, glow=glow)


def main() -> None:
    rng = np.random.default_rng(11)
    ogus = sample_ogus(rng)
    n_for_level = lambda lv, rng: (rng.integers(6, 14) if lv == 0
                                   else rng.integers(3, 6))
    for o in ogus:
        grow_children(o, rng, max_level=2, n_for_level=n_for_level)

    fig = plt.figure(figsize=(16, 9), facecolor=BG)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor(BG)
    # quantum-vacuum speckle
    ax.scatter(rng.uniform(0, W, 1100), rng.uniform(0, H, 1100),
               s=rng.uniform(0.2, 1.8, 1100), c="#9aa0c0", alpha=0.16,
               linewidths=0, zorder=0)
    drng = np.random.default_rng(11)            # reproducible colours
    for o in ogus:
        draw_node(ax, o, drng)

    # choose 'home' as deep as the foam allows: prefer a grandchild (depth 2),
    # fall back to the largest child of the largest fertile OGU.
    fertile2 = [o for o in ogus if any(c.children for c in o.children)]
    if fertile2:
        home_ogu = max(fertile2, key=lambda o: o.r)
        home_bpu = max((c for c in home_ogu.children if c.children),
                       key=lambda c: c.r)
        home = max(home_bpu.children, key=lambda c: c.r)
    else:
        home_ogu = max((o for o in ogus if o.children), key=lambda o: o.r)
        home_bpu = max(home_ogu.children, key=lambda c: c.r)
        home = home_bpu

    # magnifier inset zoomed on the home BPU (so the dot-in-dot nesting reads)
    axm = fig.add_axes([0.655, 0.04, 0.33, 0.42])
    axm.set_facecolor("#0c0c18")
    pad = home_bpu.r * 1.6
    axm.set_xlim(home_bpu.x - pad, home_bpu.x + pad)
    axm.set_ylim(home_bpu.y - pad, home_bpu.y + pad)
    axm.set_xticks([])
    axm.set_yticks([])
    for s in axm.spines.values():
        s.set_color("#39405e")
    draw_node(axm, home_ogu, np.random.default_rng(11), glow=False)
    axm.plot([home.x], [home.y], marker="o", ms=5, color="white", zorder=20)
    axm.annotate("home", (home.x, home.y), textcoords="offset points",
                 xytext=(18, 14), color="white", fontsize=12, fontweight="bold",
                 zorder=20, arrowprops=dict(arrowstyle="-", color="white", lw=1))
    axm.set_title("magnify: home, a BH inside a BHU inside an OGU",
                  color="#cdd4f0", fontsize=9, pad=3)
    con = ConnectionPatch((home_bpu.x, home_bpu.y), (home_bpu.x, home_bpu.y),
                          "data", "data", axesA=ax, axesB=axm,
                          color="white", lw=0.7, alpha=0.45, ls=(0, (4, 3)))
    fig.add_artist(con)

    fig.text(0.022, 0.07, "THE CONDENSED VOID", color="white", fontsize=23,
             fontweight="bold", alpha=0.93)
    fig.text(0.022, 0.035,
             "gravity-scaled foam  ·  log-median OGU sizes, spin-set "
             "asymmetry  ·  child ≈ 1/13 parent  ·  "
             "warm = matter, cool = antimatter  ·  rings = Cartasis membranes",
             color="#aab0d0", fontsize=10.5)

    out = os.path.join(PDF_DIR, "foam_wallpaper.pdf")
    fig.savefig(out, facecolor=BG)
    fig.savefig(out.replace(".pdf", ".png"), dpi=170, facecolor=BG)
    plt.close(fig)
    print(f"wrote {out} and .png  (OGUs={len(ogus)}, "
          f"home depth={home.level})")


if __name__ == "__main__":
    main()
