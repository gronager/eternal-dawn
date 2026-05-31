#!/usr/bin/env python3
r"""A clock's journey: in past the horizon, to the bounce, out via Hawking.

Three timescales vs black-hole mass: the clock's own (proper) infall time to the
bounce (~ M), the exterior darkening time r_s/c (~ M), and the Hawking time over
which its information returns (~ M^3). Renders figures/pdf/clock_journey.pdf.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from cartasis_sims import os_collapse as osc
from cartasis_sims import constants as k

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PDF_DIR = os.path.join(ROOT, "figures", "pdf")
OUT_DIR = os.path.join(ROOT, "sims", "output")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

YR = k.year
LANDMARKS = {
    "10 M$_\\odot$": 10 * k.M_sun,
    "10$^9$ M$_\\odot$": 1e9 * k.M_sun,
    "observable\nuniverse": 9.246e52,
}


def main() -> None:
    lines = ["A clock's journey through a black hole (proper vs exterior time)",
             "=" * 66,
             f"{'object':>22} {'infall(proper)':>16} {'darken r_s/c':>14} "
             f"{'Hawking':>14}"]
    for name, M in LANDMARKS.items():
        ti = osc.infall_proper_time(M)
        td = osc.darkening_time(M)
        th = osc.hawking_time(M)
        lines.append(f"{name.replace(chr(10),' '):>22} "
                     f"{ti:>13.2e} s {td:>11.2e} s {th/YR:>10.2e} yr")
    lines += [
        "",
        "The clock's OWN time to the bounce is finite and ~ M (0.15 ms for a",
        "stellar hole, ~22 Gyr for a universe-mass one). Its matter then seeds the",
        "baby universe; it does NOT come back as a clock. Its quantum information",
        "returns to the parent exterior smeared over the Hawking time ~ M^3 -- so",
        "the 'exit' lags the infall by a factor ~ (M^2), an astronomically larger",
        "exterior interval. There is no single 'offset': substance goes forward",
        "into the baby (~M), information comes back to the parent (~M^3).",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT_DIR, "clock_journey.txt"), "w") as f:
        f.write(text + "\n")

    M = np.logspace(30, 54, 300)
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    ax.loglog(M, [osc.infall_proper_time(m) for m in M], "C0", lw=2,
              label=r"infall (clock's proper time, $\sim\!\pi GM/c^3$)")
    ax.loglog(M, [osc.darkening_time(m) for m in M], "C2", lw=1.8, ls="--",
              label=r"exterior darkening $r_s/c$")
    ax.loglog(M, [osc.hawking_time(m) for m in M], "C3", lw=1.8, ls="-.",
              label=r"Hawking evaporation ($\sim\!M^3$)")
    ax.axhline(YR, color="0.6", lw=0.8)
    ax.text(1.2e30, YR * 1.5, "1 yr", color="0.4", fontsize=8)
    ax.axhline(13.8e9 * YR, color="0.6", lw=0.8, ls=":")
    ax.text(1.2e30, 13.8e9 * YR * 1.5, "age of universe", color="0.4", fontsize=8)
    for name, m in LANDMARKS.items():
        ax.axvline(m, color="0.85", lw=0.8, zorder=0)
        ax.text(m, 1e-7, name, rotation=90, va="bottom", ha="right", fontsize=7.5)
    ax.set_xlabel("black-hole mass  $M$  [kg]")
    ax.set_ylabel("time  [s]")
    ax.set_title("A clock's journey: proper infall ($\\sim M$) vs Hawking "
                 "return ($\\sim M^3$)", fontsize=11)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, which="major", alpha=0.2)
    fig.tight_layout()
    out_pdf = os.path.join(PDF_DIR, "clock_journey.pdf")
    fig.savefig(out_pdf)
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=130)
    plt.close(fig)
    print(f"\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
