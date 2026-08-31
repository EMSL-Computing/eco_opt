"""Generate comparison figures for all Xuezhong Wang examples."""
import os
import sys as _sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_sys.path.insert(0, os.path.dirname(__file__))

import bryson_denham_fixed_tf as ex1
import double_integrator_nonneg_control as ex2

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIGDIR, exist_ok=True)


def plot_ex1():
    sol, sys = ex1.solve()
    t_analytic = np.linspace(0, 2, 200)
    u_analytic = ex1.analytic_control(t_analytic)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(sol.time_intervals[:-1], sol.u[:, 0], ".", label="eco_opt")
    axes[0].plot(t_analytic, u_analytic, "-", label="analytic (LQ)", color="k", alpha=0.6)
    axes[0].set_xlabel("Time")
    axes[0].set_ylabel("Control (acceleration)")
    axes[0].legend()
    axes[0].grid(True)

    axes[1].plot(sol.time, sol.x[:, 0], label="position")
    axes[1].set_xlabel("Time")
    axes[1].set_ylabel("State")
    axes[1].grid(True)

    fig.suptitle(
        f"Xuezhong Wang ex.1: min-energy double integrator "
        f"(J={sol.objective:.4f}, analytic={ex1.ANALYTIC_OBJECTIVE:.4f})"
    )
    fig.tight_layout()
    out = os.path.join(FIGDIR, "xuezhong_wang_ex1.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


def plot_ex2():
    sol, sys = ex2.solve()

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].step(sol.time_intervals[:-1], sol.u[:, 0], where="post", label="eco_opt")
    axes[0].set_xlabel("Time")
    axes[0].set_ylabel("Control (acceleration)")
    axes[0].grid(True)

    axes[1].plot(sol.time, sol.x[:, 0], label="position")
    axes[1].plot(sol.time, sol.x[:, 1], label="velocity")
    axes[1].set_xlabel("Time")
    axes[1].set_ylabel("State")
    axes[1].legend()
    axes[1].grid(True)

    fig.suptitle(
        f"Xuezhong Wang ex.2: nonneg-control double integrator "
        f"(J={sol.objective:.2e}, trivial optimum J*=0)"
    )
    fig.tight_layout()
    out = os.path.join(FIGDIR, "xuezhong_wang_ex2.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    plot_ex1()
    plot_ex2()
