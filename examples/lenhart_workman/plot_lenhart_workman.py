"""Generate comparison figures for all Lenhart & Workman examples."""
import os
import sys as _sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_sys.path.insert(0, os.path.dirname(__file__))

import linear_bangbang_ex17_2 as ex172
import lq_tracking_ex17_3 as ex173
import min_time_trolley_ex20_4 as ex204

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIGDIR, exist_ok=True)


def plot_ex172():
    sol, sys = ex172.solve()
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].step(sol.time_intervals[:-1], sol.u[:, 0], where="post", label="eco_opt")
    ax[0].axvline(ex172.ANALYTIC_SWITCH_TIME, color="k", ls="--", label="analytic switch")
    ax[0].set_xlabel("Time")
    ax[0].set_ylabel("Control u")
    ax[0].legend()
    ax[0].grid(True)
    ax[1].plot(sol.time, sol.x[:, 0])
    ax[1].set_xlabel("Time")
    ax[1].set_ylabel("State x")
    ax[1].grid(True)
    fig.suptitle(f"Lenhart-Workman ex.17.2: bang-bang (J={sol.objective:.4f})")
    fig.tight_layout()
    out = os.path.join(FIGDIR, "lenhart_workman_ex17_2.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


def plot_ex173():
    sol, sys = ex173.solve()
    t_ref = np.linspace(0, 2, 200)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(sol.time_intervals[:-1], sol.u[:, 0], ".-", label="eco_opt")
    ax[0].set_xlabel("Time")
    ax[0].set_ylabel("Control u")
    ax[0].grid(True)
    ax[1].plot(sol.time, sol.x[:, 0], label="x(t)")
    ax[1].plot(t_ref, t_ref ** 2, "k--", alpha=0.6, label="reference t^2")
    ax[1].set_xlabel("Time")
    ax[1].set_ylabel("State")
    ax[1].legend()
    ax[1].grid(True)
    fig.suptitle(f"Lenhart-Workman ex.17.3: LQ tracking (J={sol.objective:.4f})")
    fig.tight_layout()
    out = os.path.join(FIGDIR, "lenhart_workman_ex17_3.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


def plot_ex204():
    sol, sys = ex204.solve()
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    ax[0].plot(sol.time, sol.x[:, 0])
    ax[0].set_xlabel("Time")
    ax[0].set_ylabel("Position")
    ax[0].grid(True)
    ax[1].plot(sol.time, sol.x[:, 1])
    ax[1].set_xlabel("Time")
    ax[1].set_ylabel("Velocity")
    ax[1].grid(True)
    ax[2].step(sol.time_intervals[:-1], sol.u[:, 0], where="post")
    ax[2].set_xlabel("Time")
    ax[2].set_ylabel("Acceleration (control)")
    ax[2].grid(True)
    fig.suptitle(
        f"Lenhart-Workman ex.20.4: min-time stop "
        f"(tf={sol.tf:.4f}, analytic={ex204.analytic_tf():.4f})"
    )
    fig.tight_layout()
    out = os.path.join(FIGDIR, "lenhart_workman_ex20_4.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    plot_ex172()
    plot_ex173()
    plot_ex204()
