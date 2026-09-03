"""Cycloid overlay for the brachistochrone example."""
import os
import sys as _sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_sys.path.insert(0, os.path.dirname(__file__))
import brachistochrone as br

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIGDIR, exist_ok=True)


def plot_brachistochrone():
    sol, sys = br.solve()
    t_fine = np.linspace(0.0, br.ANALYTIC_TF, 400)
    xa, ya, va, tha = br.analytic_trajectories(t_fine)

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.2))

    ax = axes[0, 0]
    ax.plot(sol.x[:, 0], sol.x[:, 1], "o", ms=3, color="C0", label="eco_opt")
    ax.plot(xa, ya, "-", color="k", alpha=0.55, label="cycloid")
    ax.plot([0, br.X_F], [0, br.Y_F], "k+", ms=10)
    ax.set_xlabel("x (horizontal)")
    ax.set_ylabel("y (downward)")
    ax.invert_yaxis()
    ax.set_aspect("equal", adjustable="box")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax.set_title("wire shape")

    ax = axes[0, 1]
    ax.step(sol.time_intervals[:-1], sol.u[:, 0], where="post", color="C2", label="eco_opt $\\theta$")
    ax.plot(t_fine, tha, "-", color="k", alpha=0.55, label="analytic $\\theta=\\varphi/2$")
    ax.axhline(0.5 * np.pi, color="0.5", ls=":", lw=1)
    ax.set_xlabel("time")
    ax.set_ylabel("heading $\\theta$ from vertical")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1, 0]
    ax.plot(sol.time, sol.x[:, 2], "o", ms=3, color="C3", label="eco_opt $v$")
    ax.plot(t_fine, va, "-", color="k", alpha=0.55, label="analytic $v$")
    ax.set_xlabel("time")
    ax.set_ylabel("speed")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    ax.plot(sol.time, sol.x[:, 0], "o", ms=3, color="C0", label="eco_opt $x$")
    ax.plot(sol.time, sol.x[:, 1], "s", ms=3, color="C1", label="eco_opt $y$")
    ax.plot(t_fine, xa, "-", color="C0", alpha=0.45)
    ax.plot(t_fine, ya, "-", color="C1", alpha=0.45)
    ax.set_xlabel("time")
    ax.set_ylabel("position")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.suptitle(
        "Brachistochrone (Bernoulli): min-time cycloid\n"
        f"$t_f$={sol.tf:.5f} (analytic $\\pi$={br.ANALYTIC_TF:.5f})",
        fontsize=11,
    )
    fig.tight_layout()
    out = os.path.join(FIGDIR, "brachistochrone.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    plot_brachistochrone()
