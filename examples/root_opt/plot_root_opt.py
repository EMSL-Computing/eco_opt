"""Compare three scenarios (base, low mineral N, high SOM C:N) over a
120-day growing season.
"""
import copy
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

SCENARIOS = {}


def run_scenarios(control_intervals=50, polynomial_degree=3):
    results = {}

    # base
    results["base"] = m.solve(control_intervals=control_intervals, polynomial_degree=polynomial_degree)

    # low mineral N
    state_low_n = copy.deepcopy(m.STATE_INIT)
    state_low_n["inorgN"] = 10.0
    results["low min. N"] = m.solve(
        state_init=state_low_n, control_intervals=control_intervals, polynomial_degree=polynomial_degree
    )

    # high SOM C:N (SOMN reduced -> higher SOM C:N ratio)
    state_high_cn = copy.deepcopy(m.STATE_INIT)
    state_high_cn["SOMN"] = state_high_cn["SOMC"] / 20.0
    results["high SOM C:N"] = m.solve(
        state_init=state_high_cn, control_intervals=control_intervals, polynomial_degree=polynomial_degree
    )

    return results


def plot_scenarios(results):
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    colors = {"base": "tab:gray", "low min. N": "tab:orange", "high SOM C:N": "tab:green"}

    for name, (sol, sys) in results.items():
        c = colors[name]
        axes[0, 0].plot(sol.time, sol.x[:, 0], color=c, label=name)  # leafC
        axes[0, 1].plot(sol.time, sol.x[:, 1], color=c, label=name)  # rootC
        axes[0, 2].plot(sol.time, sol.x[:, 8], color=c, label=name)  # inorgN
        axes[1, 0].plot(sol.time_intervals[:-1], sol.u[:, 0], color=c, label=name)  # rootCSupply control
        axes[1, 1].step(sol.time_intervals[:-1], sol.u[:, 1], where="post", color=c, label=name)  # exu_fraction control
        axes[1, 2].bar(name, sol.objective, color=c)

    axes[0, 0].set_title("Leaf C [gC/m2]")
    axes[0, 1].set_title("Root C [gC/m2]")
    axes[0, 2].set_title("Mineral N [gN/m2]")
    axes[1, 0].set_title("Control: root C supply [gC/m2/d]")
    axes[1, 1].set_title("Control: exudation fraction [-]")
    axes[1, 2].set_title("Objective (cumulative leaf growth) [gC/m2]")

    for ax in axes[0]:
        ax.set_xlabel("time [d]")
        ax.grid(True)
    for ax in axes[1, :2]:
        ax.set_xlabel("time [d]")
        ax.grid(True)
    axes[0, 0].legend(fontsize=8)

    fig.suptitle("root_opt: plant C-N allocation, 3-scenario comparison (eco_opt port)")
    fig.tight_layout()
    out = os.path.join(FIGDIR, "root_opt_scenarios.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    results = run_scenarios()
    for name, (sol, sys) in results.items():
        print(f"{name}: success={sol.success} objective={sol.objective:.4f}")
    plot_scenarios(results)
