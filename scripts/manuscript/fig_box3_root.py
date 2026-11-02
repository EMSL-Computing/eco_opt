"""Coordinated allocation in the published plant--soil model.

(a) Coupled processes and the two controls, with the growth objective
separate from the schematic. (b) Belowground C supply relative to net
assimilation. (c) Exudation fraction. (d) Leaf C accumulation.

Model and scenarios adapted from Chakrawal, Mooney and Colombi (2026).
The numerical formulation is unchanged; this figure uses two optimal
solutions and does not compute shadow-price or perturbation diagnostics.
"""
from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path

import casadi as cs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import (
    apply_style, panel_label, save, shared_legend, TEXT_W,
    STATE_COLOR, CONTROL_COLOR, PAYOFF_COLOR, RULE_LS,
)
from models import plant_soil as ro

apply_style()
P = ro.PARAMS
N_HIGH, N_LOW = 20.0, 10.0
N_INTERVALS = 100
LS = {N_HIGH: "-", N_LOW: (0, (1.2, 1.2))}
LW = 1.5
BOUND_COLOR = "0.6"
_x = cs.MX.sym("x", 9)
_u = cs.MX.sym("u", 2)
_, _y = ro.plant_soil_microbe_nested(0, _x, _u)
NET_ASSIMILATION = cs.Function("net_assimilation", [_x, _u], [_y.Net_leafC_assimilation])


def solve(inorgN):
    state_init = copy.deepcopy(ro.STATE_INIT)
    state_init["inorgN"] = inorgN
    sol, _ = ro.solve(state_init=state_init, control_intervals=N_INTERVALS)
    assert sol.success, sol.return_status
    return sol


def analyse(sol):
    tu = sol.time_intervals[:-1]
    states = np.column_stack([np.interp(tu, sol.time, sol.x[:, i]) for i in range(9)])
    assimilation = np.array([float(NET_ASSIMILATION(x, u))
                             for x, u in zip(states, sol.u)])
    return {"t": sol.time, "tu": tu, "S_frac": sol.u[:, 0] / assimilation,
            "f": sol.u[:, 1], "leafC": sol.x[:, 0], "J": float(sol.objective)}


def draw_schematic(ax):
    """Arrange the resource pathways without crossing arrows or overlaid text."""
    ax.set_xlim(-0.03, 1.07)
    ax.set_ylim(-0.30, 1.03)
    ax.axis("off")
    grey = "0.35"

    def pool(x, y, w, h, text, color):
        ax.add_patch(FancyBboxPatch((x, y), w, h, fc="white", ec=color, lw=1.0,
                                    boxstyle="round,pad=0.0,rounding_size=0.025"))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=7.5)

    def arrow(p, q, color, lw=1.0):
        ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=8,
                                     color=color, lw=lw, shrinkA=1, shrinkB=1))

    def label(x, y, text, color=grey, **kw):
        ax.text(x, y, text, fontsize=7, color=color, **kw)

    # Plant pools form a single column; soil processes occupy the right.
    pool(0.04, 0.77, 0.22, 0.19, "leaves", PAYOFF_COLOR)
    pool(0.04, 0.43, 0.22, 0.19, "roots", grey)
    pool(0.04, 0.09, 0.22, 0.19, "mineral N", STATE_COLOR)
    pool(0.50, 0.43, 0.20, 0.19, "exudates", grey)
    pool(0.80, 0.77, 0.22, 0.19, "soil organic C, N", STATE_COLOR)
    pool(0.80, 0.09, 0.22, 0.19, "microbes", grey)

    # Two decisions, with C investment and returning N shown side by side.
    arrow((0.105, 0.77), (0.105, 0.62), CONTROL_COLOR, lw=1.5)
    label(0.085, 0.695, r"C supply $S$", CONTROL_COLOR, ha="right", va="center")
    arrow((0.195, 0.62), (0.195, 0.77), STATE_COLOR)
    label(0.215, 0.695, "N for assimilation", STATE_COLOR, ha="left", va="center")
    arrow((0.26, 0.525), (0.50, 0.525), CONTROL_COLOR, lw=1.5)
    label(0.38, 0.55, r"exudation $f e_RS$", CONTROL_COLOR, ha="center", va="bottom")

    # Acquisition through roots and microbial processing close the feedback.
    arrow((0.15, 0.28), (0.15, 0.43), STATE_COLOR)
    label(0.17, 0.355, "root N uptake", STATE_COLOR, ha="left", va="center")
    arrow((0.70, 0.46), (0.84, 0.28), grey)
    label(0.76, 0.49, "C uptake", ha="center", va="bottom")
    arrow((0.91, 0.77), (0.91, 0.28), grey)
    label(0.93, 0.525, "C, N uptake", ha="left", va="center")
    arrow((0.80, 0.185), (0.26, 0.185), STATE_COLOR)
    label(0.53, 0.225, "net N mineralization", STATE_COLOR, ha="center", va="bottom")

    # The objective has its own strip below the complete process diagram.
    ax.plot([0.04, 1.02], [-0.035, -0.035], color="0.82", lw=0.6)
    ax.text(0.53, -0.19, r"Objective   maximize cumulative leaf growth $J=\int_0^{120}G_L\,dt$",
            ha="center", va="center", fontsize=7.5, color=PAYOFF_COLOR)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot-data", type=Path,
                        help="Save/reuse numerical trajectories while adjusting figure layout.")
    args = parser.parse_args()
    if args.plot_data is not None and args.plot_data.exists():
        with np.load(args.plot_data) as data:
            recs = {N0: {key: data[f"{N0}_{key}"] for key in
                          ("t", "tu", "S_frac", "f", "leafC", "J")}
                    for N0 in (N_HIGH, N_LOW)}
    else:
        recs = {N0: analyse(solve(N0)) for N0 in (N_HIGH, N_LOW)}
        if args.plot_data is not None:
            np.savez(args.plot_data, **{f"{N0}_{key}": value
                                      for N0, rec in recs.items()
                                      for key, value in rec.items()})
    fig = plt.figure(figsize=(TEXT_W, 3.35))
    gs = fig.add_gridspec(2, 3, wspace=0.68, hspace=0.35, height_ratios=[1.35, 1])
    ax = {"a": fig.add_subplot(gs[0, :])}
    ax.update({k: fig.add_subplot(gs[1, j]) for j, k in enumerate("bcd")})
    draw_schematic(ax["a"])
    panel_label(ax["a"], "a", x=0.02, y=0.98)

    a_ = ax["b"]
    for N0, rec in recs.items():
        a_.step(rec["tu"], rec["S_frac"], where="post", color=CONTROL_COLOR, ls=LS[N0], lw=LW)
    a_.set(xlim=(0, 120), ylim=(0, 1), xlabel="time (d)",
           ylabel="belowground C supply\n" + r"$S/A_{\mathrm{net}}$")
    panel_label(a_, "b")

    a_ = ax["c"]
    for bound in (P["min_exudation"], P["max_exudation"]):
        a_.axhline(bound, color=BOUND_COLOR, ls=RULE_LS, lw=0.8)
    for N0, rec in recs.items():
        a_.step(rec["tu"], rec["f"], where="post", color=CONTROL_COLOR, ls=LS[N0], lw=LW)
    a_.set(xlim=(0, 120), ylim=(0, 0.62), xlabel="time (d)",
           ylabel="exudation fraction $f$")
    panel_label(a_, "c")

    a_ = ax["d"]
    for N0, rec in recs.items():
        a_.plot(rec["t"], rec["leafC"], color=PAYOFF_COLOR, ls=LS[N0], lw=LW)
    a_.set(xlim=(0, 120), ylim=(0, None), xlabel="time (d)",
           ylabel=r"leaf C (g C m$^{-2}$)")
    panel_label(a_, "d")

    shared_legend(fig, [
        (rf"initial N: {N_HIGH:.0f} g N m$^{{-2}}$", dict(color="0.2", lw=LW, ls=LS[N_HIGH])),
        (rf"initial N: {N_LOW:.0f} g N m$^{{-2}}$", dict(color="0.2", lw=LW, ls=LS[N_LOW])),
        ("bounds on $f$", dict(color=BOUND_COLOR, lw=0.8, ls=RULE_LS)),
    ], y=0.95)
    save(fig, "fig_box3_root")
    for N0, rec in recs.items():
        print(f"N0={N0:.0f}: J*={rec['J']:.2f} g C m-2; leaf C(T)={rec['leafC'][-1]:.2f}")


if __name__ == "__main__":
    main()
