"""Generate comparison figures for the Goddard rocket free-final-time problem.

Run:  python examples/goddard_benchmark/plot_goddard_benchmark.py
Output: figs/goddard_states_comparison.png, figs/goddard_control_comparison.png
"""
import os
import sys as _sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.io import loadmat

HERE = os.path.dirname(os.path.abspath(__file__))
_sys.path.insert(0, HERE)
import goddard_benchmark as gb

ROOT = os.path.dirname(os.path.dirname(HERE))
REF_DIR = os.path.join(ROOT, "work", "goddard_benchmark")
FIG_DIR = os.path.join(HERE, "figs")


def main():
    os.makedirs(FIG_DIR, exist_ok=True)

    sol, sys = gb.solve(control_intervals=59, polynomial_degree=5, tf_guess=200.0)
    assert sol.success, sol.return_status

    yop = loadmat(os.path.join(REF_DIR, "yop_goddard.mat"))
    yop_t, yop_x, yop_u = yop["yop_t"].ravel(), yop["yop_x"], yop["yop_u"].ravel()

    casadi_path = os.path.join(REF_DIR, "casadi_goddard.mat")
    casadi_ref = loadmat(casadi_path) if os.path.exists(casadi_path) else None

    propt_path = os.path.join(REF_DIR, "propt_goddard.mat")
    propt_ref = loadmat(propt_path) if os.path.exists(propt_path) else None

    # -- States figure --------------------------------------------------
    fig, axes = plt.subplots(3, 1, figsize=(7, 9), sharex=True)
    labels = ["velocity [m/s]", "height [m]", "mass [kg]"]
    for i, ax in enumerate(axes):
        ax.plot(sol.time, sol.x[:, i], "o-", color="C0", ms=3, label="eco_opt")
        ax.plot(yop_t, yop_x[:, i], "x--", color="C1", ms=5, label="reference")
        if casadi_ref is not None:
            try:
                ax.plot(casadi_ref["casadi_t"].ravel(), casadi_ref["casadi_x"][:, i],
                        ".:", color="C2", label="CasADi (reference)")
            except KeyError:
                pass
        if propt_ref is not None:
            try:
                ax.plot(propt_ref["propt_t"].ravel(), propt_ref["propt_x"][:, i],
                        ".:", color="C3", label="PROPT (reference)")
            except KeyError:
                pass
        ax.set_ylabel(labels[i])
        ax.grid(alpha=0.3)
    axes[0].legend(loc="best", fontsize=8)
    axes[-1].set_xlabel("time [s]")
    fig.suptitle("Goddard rocket (free tf): eco_opt vs. reference states")
    fig.tight_layout()
    states_path = os.path.join(FIG_DIR, "goddard_states_comparison.png")
    fig.savefig(states_path, dpi=150)
    plt.close(fig)

    # -- Control figure ---------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.step(sol.time_intervals, np.append(sol.u[:, 0], sol.u[-1, 0]), where="post",
            color="C0", label="eco_opt")
    ax.step(yop_t, yop_u, where="post", color="C1", linestyle="--", label="reference")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("fuel mass flow F (control)")
    ax.set_title("Goddard rocket (free tf): eco_opt vs. reference control")
    ax.legend(loc="best", fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    control_path = os.path.join(FIG_DIR, "goddard_control_comparison.png")
    fig.savefig(control_path, dpi=150)
    plt.close(fig)

    print(f"objective (final height): eco_opt={sol.objective:.6f}, "
          f"yop={yop_x[-1,1]:.6f}")
    print(f"final time: eco_opt={sol.tf:.6f}, yop={yop_t[-1]:.6f}")
    print("Wrote:", states_path)
    print("Wrote:", control_path)


if __name__ == "__main__":
    main()
