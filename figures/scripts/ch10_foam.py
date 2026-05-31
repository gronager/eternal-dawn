#!/usr/bin/env python3
r"""Foam wallpaper (grayscale): the gravity-scaled condensed void.

A self-consistent illustration of the thesis (Chapters 5, 10), now in grayscale
so the structure -- not colour -- carries the meaning, and print-safe:

* Gray void. Every universe is a gray disc on a gray background; the Cartasis
  membranes (peripheries) carry all the contrast.
* OG universes (OGUs == BHU_0) get a bold DOUBLE periphery -- a black ring with a
  white ring inside it -- so they read as the top-level condensations.
* Nested universes get a single periphery, BLACK or WHITE by chirality
  (black = matter, white = antimatter). The bias is set per bounce (CPT-even),
  random sign -- not inherited.
* Sizes: OGU radii are log-normal (defined median): births outrun runaway growth,
  so many small + few large. Nesting is self-similar at diameter ratio ~1/13
  (gravity-scaled, NOT the ~1e-12 mass ratio).
* 'home' is marked deep in the foam. The framework's own result (Ch. 5) is that
  we are NOT an OGU: dark matter + dark energy require a parent, so our depth is
  n >= 1 -- a range, to be tightened by the DM/DE amplitudes and C. The drawn
  depth is illustrative; the label says n>=1.

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

# grayscale palette
VOID = "0.50"          # mid-gray background void
FILL = "0.50"          # universe interiors: gray (membranes carry contrast)
BLACK = "0.05"
WHITE = "0.97"
NEST = 1.0 / 13.0      # child / parent diameter ratio (log-gravity scale)
W, H = 1.78, 1.0       # 16:9 void canvas


@dataclass
class Node:
    x: float
    y: float
    r: float
    level: int
    matter: bool
    children: list = field(default_factory=list)


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
    radii = np.sort(rng.lognormal(mean=np.log(median_r), sigma=sigma,
                                   size=n * 3))[::-1]
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


def draw_node(ax, node, rng):
    """Gray disc; membrane drawn per the grayscale spec."""
    z = node.level + 1
    # gray interior (very slight level tint so deep nests don't vanish)
    fill = str(min(0.62, 0.50 + 0.04 * node.level))
    if node.level == 0:
        # OGU: gray fill + DOUBLE periphery (black ring, white ring inside)
        ax.add_patch(Circle((node.x, node.y), node.r, facecolor=fill,
                            edgecolor=BLACK, lw=1.7, zorder=z))
        ax.add_patch(Circle((node.x, node.y), node.r * 0.965, facecolor="none",
                            edgecolor=WHITE, lw=0.9, zorder=z + 0.5))
        # spin-axis tick on larger OGUs (median viable spin -> sets C)
        if node.r > 0.045:
            ang = rng.random() * np.pi
            dx, dy = node.r * 0.7 * np.cos(ang), node.r * 0.7 * np.sin(ang)
            ax.plot([node.x - dx, node.x + dx], [node.y - dy, node.y + dy],
                    color=BLACK, lw=0.5, alpha=0.45, zorder=z + 1)
    else:
        # nested: single periphery, black=matter / white=antimatter
        edge = BLACK if node.matter else WHITE
        lw = max(0.4, 1.1 - 0.35 * (node.level - 1))
        ax.add_patch(Circle((node.x, node.y), node.r, facecolor=fill,
                            edgecolor=edge, lw=lw, zorder=z))
    for ch in node.children:
        draw_node(ax, ch, rng)


def main() -> None:
    rng = np.random.default_rng(11)
    ogus = sample_ogus(rng)
    n_for_level = lambda lv, rng: (rng.integers(6, 14) if lv == 0
                                   else rng.integers(3, 6))
    for o in ogus:
        grow_children(o, rng, max_level=2, n_for_level=n_for_level)

    fig = plt.figure(figsize=(16, 9), facecolor=VOID)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor(VOID)
    # quantum-vacuum speckle: faint darker + lighter gray flecks
    sp = rng.uniform(0, W, 1300), rng.uniform(0, H, 1300)
    shade = rng.choice([0.40, 0.60], size=1300)
    ax.scatter(*sp, s=rng.uniform(0.2, 1.6, 1300),
               c=[str(s) for s in shade], alpha=0.30, linewidths=0, zorder=0)
    drng = np.random.default_rng(11)
    for o in ogus:
        draw_node(ax, o, drng)

    # 'home' as deep as the foam allows (illustrative; true depth n>=1, a range)
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

    # magnifier inset zoomed on the home BPU
    axm = fig.add_axes([0.655, 0.04, 0.33, 0.42])
    axm.set_facecolor("0.44")
    pad = home_bpu.r * 1.6
    axm.set_xlim(home_bpu.x - pad, home_bpu.x + pad)
    axm.set_ylim(home_bpu.y - pad, home_bpu.y + pad)
    axm.set_xticks([])
    axm.set_yticks([])
    for s in axm.spines.values():
        s.set_color(BLACK)
        s.set_linewidth(1.0)
    draw_node(axm, home_ogu, np.random.default_rng(11))
    axm.plot([home.x], [home.y], marker="o", ms=5, color=BLACK,
             markeredgecolor=WHITE, markeredgewidth=0.8, zorder=20)
    axm.annotate("home  ($n\\geq1$: not an OGU)", (home.x, home.y),
                 textcoords="offset points", xytext=(16, 14), color=BLACK,
                 fontsize=11, fontweight="bold", zorder=20,
                 arrowprops=dict(arrowstyle="-", color=BLACK, lw=1))
    axm.set_title("magnify: home, a BH-universe nested inside an OGU",
                  color=BLACK, fontsize=9, pad=3)
    con = ConnectionPatch((home_bpu.x, home_bpu.y), (home_bpu.x, home_bpu.y),
                          "data", "data", axesA=ax, axesB=axm,
                          color=BLACK, lw=0.8, alpha=0.5, ls=(0, (4, 3)))
    fig.add_artist(con)

    fig.text(0.022, 0.07, "THE CONDENSED VOID", color=BLACK, fontsize=23,
             fontweight="bold", alpha=0.92)
    fig.text(0.022, 0.035,
             "gravity-scaled foam  ·  log-median OGU sizes  ·  child ≈ 1/13 "
             "parent  ·  bold double ring = OGU membrane  ·  "
             "black = matter, white = antimatter membranes",
             color=BLACK, fontsize=10.5, alpha=0.82)

    out = os.path.join(PDF_DIR, "foam_wallpaper.pdf")
    fig.savefig(out, facecolor=VOID)
    fig.savefig(out.replace(".pdf", ".png"), dpi=170, facecolor=VOID)
    plt.close(fig)
    print(f"wrote {out} and .png  (OGUs={len(ogus)}, home depth={home.level})")


if __name__ == "__main__":
    main()
