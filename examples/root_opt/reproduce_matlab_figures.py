"""Reproduce (as closely as matplotlib allows) the two MATLAB figures that
are actually exported by ``external/root_opt/Figure2_and_3.m``:

  * Figure2_revised.png -- 5x3 tiled layout, panels (A)-(J): Anet time
    series + cumulative bar, then leaf/root/exudation growth-rate RATIOS
    (G_L/A_net, S/A_net, G_R/S_net, E/S_net) time series + cumulative bars,
    for the 3 scenarios (base / low min. N / high SOM C:N).
  * Figure3.png -- 4x1 tiled layout, panels (A)-(D): root N uptake rate,
    root growth rate, root exudation rate, net N mineralisation (phiN)
    time series for the 3 scenarios.

Matches (as closely as matplotlib/GridSpec allow): panel layout and
span, panel order, axis labels (incl. italic/subscript via mathtext),
copper(3) 3-color scenario palette, line widths, grid/box on every
axes, axis limits, (A)-(J)/(A)-(D) panel titles top-left, legend on the
first axes, and hidden x-tick labels on non-bottom rows -- i.e. the
same *figure*, not just the same data, as ``Figure2_and_3.m`` produces.
"""
import os
import sys as _sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_sys.path.insert(0, os.path.dirname(__file__))
import plant_soil_nitrogen_allocation as m

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIGDIR, exist_ok=True)

SCENARIOS = ["base", "low min. N", "high SOM C:N"]
# matplotlib ships the same 'copper' colormap MATLAB's copper(3) uses.
COLORS = matplotlib.colormaps["copper"](np.linspace(0, 1, 3))


def run_scenarios(terminal_time=120.0, control_intervals=100, polynomial_degree=3):
    """Solve the 3 scenarios (matches MATLAB's control_intervals=100)."""
    results = {}
    results["base"] = m.solve(terminal_time=terminal_time,
                               control_intervals=control_intervals,
                               polynomial_degree=polynomial_degree)

    state_low_n = dict(m.STATE_INIT)
    state_low_n["inorgN"] = 7.5
    results["low min. N"] = m.solve(terminal_time=terminal_time, state_init=state_low_n,
                                     control_intervals=control_intervals,
                                     polynomial_degree=polynomial_degree)

    state_high_cn = dict(m.STATE_INIT)
    state_high_cn["SOMN"] = state_high_cn["SOMC"] / 20.0
    results["high SOM C:N"] = m.solve(terminal_time=terminal_time, state_init=state_high_cn,
                                       control_intervals=control_intervals,
                                       polynomial_degree=polynomial_degree)
    return results


def _signals(sol, sys):
    """Sample every derived variable used by the two figures on the coarse
    (interval-boundary) time grid, matching MATLAB's per-row-of-df access."""
    y = sys.y
    out = {"time": sol.time}
    for name in ["Anet", "leaf_growth_rate", "rootCSupply", "root_growth_rate",
                 "root_exu", "Net_rootC_assimilation", "root_N_uptake", "phiN"]:
        out[name] = sol.signal(getattr(y, name))
    return out


def make_figure2_revised(data):
    fig = plt.figure(figsize=(1100 / 100, 1100 / 100), dpi=100)
    fig.patch.set_facecolor("w")
    gs = fig.add_gridspec(5, 3, hspace=0.45, wspace=0.45)

    axA = fig.add_subplot(gs[0, 0:2]); axB = fig.add_subplot(gs[0, 2])
    axC = fig.add_subplot(gs[1, 0:2]); axD = fig.add_subplot(gs[1, 2])
    axE = fig.add_subplot(gs[2, 0:2]); axF = fig.add_subplot(gs[2, 2])
    axG = fig.add_subplot(gs[3, 0:2]); axH = fig.add_subplot(gs[3, 2])
    axI = fig.add_subplot(gs[4, 0:2]); axJ = fig.add_subplot(gs[4, 2])
    ax = [axA, axB, axC, axD, axE, axF, axG, axH, axI, axJ]
    for a in ax:
        a.grid(True)

    alpha_val = 0.95

    # ---- Row 1: Anet time series + cumulative bar --------------------
    cumA = []
    for i, name in enumerate(SCENARIOS):
        df = data[name]
        axA.plot(df["time"], df["Anet"], lw=2, color=COLORS[i], label=name)
        cumA.append(np.trapz(df["Anet"], df["time"]))
    axB.bar(SCENARIOS, cumA, color=COLORS, alpha=alpha_val, edgecolor=COLORS)
    axA.set_ylabel(r"$A_{net}$ [gC m$^{-2}$ d$^{-1}$]")
    axB.set_ylabel(r"$\int A_{net}$ [gC m$^{-2}$]")

    # ---- Rows 2-5: ratio time series + cumulative bars ----------------
    yvar_ratio = ["leaf_growth_rate", "rootCSupply", "root_growth_rate", "root_exu"]
    ratio_den = ["Anet", "Anet", "Net_rootC_assimilation", "Net_rootC_assimilation"]
    ylabels_ratio = [r"$G_{L}/A_{net}$", r"$S/A_{net}$", r"$G_{R}/S_{net}$", r"$E/S_{net}$"]
    ylims_ratio = [1, 1, 1, 0.6]
    ax_ratio = [axC, axE, axG, axI]
    ax_bar = [axD, axF, axH, axJ]

    for j in range(4):
        temp = []
        for i, name in enumerate(SCENARIOS):
            df = data[name]
            temp.append(np.trapz(df[yvar_ratio[j]], df["time"]))
        ax_bar[j].bar(SCENARIOS, temp, color=COLORS, alpha=alpha_val, edgecolor=COLORS)
        base_label = ylabels_ratio[j].split("/")[0].lstrip("$")
        ax_bar[j].set_ylabel(r"$\int$ " + base_label + " [gC m$^{-2}$]")

        for i, name in enumerate(SCENARIOS):
            df = data[name]
            line_color = COLORS[i] * alpha_val + (1 - alpha_val) * np.array([1, 1, 1, 1])
            ax_ratio[j].plot(df["time"], np.array(df[yvar_ratio[j]]) / np.array(df[ratio_den[j]]),
                              lw=2, color=line_color, label=name)
        ax_ratio[j].set_ylabel(ylabels_ratio[j] + " [-]")
        ax_ratio[j].set_ylim(0, ylims_ratio[j])

    for a in [axC, axE, axG, axI]:
        a.axhline(0.2 if a is not axE else 0.8, ls="--", color=[0.25, 0.25, 0.25], lw=1)
    axI.axhline(m.PARAMS["min_exudation"], ls="--", color=[0.25, 0.25, 0.25], lw=1)
    axI.axhline(m.PARAMS["max_exudation"], ls="--", color=[0.25, 0.25, 0.25], lw=1)

    axE.set_ylim(0, 1)
    for a in [axD, axF]:
        a.set_ylim(0, 900)
    for a in [axH, axJ]:
        a.set_ylim(0, 300)

    for a in [axA, axB, axC, axD, axE, axF, axG, axH]:
        a.set_xticklabels([])
    axI.set_xlabel("time [d]")
    axJ.set_xlabel("")

    leg = axA.legend(loc="upper left", fontsize=10, ncol=2, frameon=False)

    titles = "ABCDEFGHIJ"
    for k, a in enumerate(ax):
        a.set_title(f"({titles[k]})", loc="left", fontweight="normal", fontsize=14,
                    x=0.0, y=1.0)
        a.tick_params(labelsize=10)
        for spine in a.spines.values():
            spine.set_linewidth(0.5)

    out = os.path.join(FIGDIR, "Figure2_revised.png")
    fig.savefig(out, dpi=150, facecolor="w")
    plt.close(fig)
    print("wrote", out)


def make_figure3(data):
    fig, ax = plt.subplots(4, 1, figsize=(700 / 100, 850 / 100), dpi=100, sharex=False)
    fig.patch.set_facecolor("w")
    ax2, ax3, ax4, ax5 = ax
    for a in ax:
        a.grid(True)

    for i, name in enumerate(SCENARIOS):
        df = data[name]
        ax2.plot(df["time"], df["root_N_uptake"], lw=2, color=COLORS[i], label=name)
        ax3.plot(df["time"], df["root_growth_rate"], lw=2, color=COLORS[i], label=name)
        ax4.plot(df["time"], df["root_exu"], lw=2, color=COLORS[i], label=name)
        ax5.plot(df["time"], df["phiN"], lw=2, color=COLORS[i], label=name)

    ax3.set_ylabel("root growth rate\n" + r"$G_R$ [gC m$^{-2}$ d$^{-1}$]")
    ax4.set_ylabel("exudation rate\n" + r"$E$ [gC m$^{-2}$ d$^{-1}$]")
    ax2.set_ylabel("root N uptake rate\n" + r"$U_N$ [gN m$^{-2}$ d$^{-1}$]")
    ax5.set_ylabel("net mineralisation\n" + r"$\phi_N$ [gN m$^{-2}$ d$^{-1}$]")
    ax5.set_xlabel("time [d]")

    for a in [ax2, ax3, ax4]:
        a.set_xticklabels([])

    titles = "ABCD"
    for k, a in enumerate(ax):
        a.set_title(f"({titles[k]})", loc="left", fontweight="normal", fontsize=14, x=0.0, y=1.0)
        a.tick_params(labelsize=11)
        for spine in a.spines.values():
            spine.set_linewidth(0.5)

    leg = ax2.legend(loc="best", fontsize=13, ncol=3, frameon=True)

    fig.tight_layout()
    out = os.path.join(FIGDIR, "Figure3.png")
    fig.savefig(out, dpi=150, facecolor="w")
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    results = run_scenarios()
    data = {}
    for name, (sol, sys) in results.items():
        assert sol.success, f"{name} scenario failed to converge"
        data[name] = _signals(sol, sys)
        print(f"{name}: objective={sol.objective:.4f}")

    make_figure2_revised(data)
    make_figure3(data)
