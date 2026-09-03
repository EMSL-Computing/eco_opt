"""Figures for linear and quadratic QSS examples (ODE vs DAE)."""
import os
import sys as _sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_sys.path.insert(0, os.path.dirname(__file__))
import linear_qss as lin
import quadratic_qss as quad

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIGDIR, exist_ok=True)


def _plot_pair(name, sol_ode, sol_dae, M_of):
    fig, axes = plt.subplots(3, 1, figsize=(7.0, 7.4), sharex=True)
    axes[0].plot(sol_ode.time, sol_ode.x[:, 0], "-", color="C0", label="ODE (reduced)")
    axes[0].plot(sol_dae.time, sol_dae.x[:, 0], "o", ms=3, color="C1", label="DAE")
    axes[0].set_ylabel("$C$")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3)

    M_ode = np.array([float(M_of(c, u)) for c, u in zip(sol_ode.x[:, 0], np.append(sol_ode.u[:, 0], sol_ode.u[-1, 0]))])
    axes[1].plot(sol_ode.time, M_ode, "-", color="C0")
    axes[1].plot(sol_dae.time, sol_dae.z[:, 0], "o", ms=3, color="C1")
    axes[1].set_ylabel("$M$")
    axes[1].grid(alpha=0.3)

    t_u = sol_ode.time_intervals
    axes[2].step(t_u[:-1], sol_ode.u[:, 0], where="post", color="C0")
    axes[2].step(sol_dae.time_intervals[:-1], sol_dae.u[:, 0], where="post", color="C1", ls="--")
    axes[2].set_ylabel("$u$")
    axes[2].set_xlabel("time")
    axes[2].grid(alpha=0.3)

    fig.suptitle(
        f"{name}: reduced ODE vs DAE\n"
        f"J_ode={sol_ode.objective:.5f}, J_dae={sol_dae.objective:.5f}",
        fontsize=11,
    )
    fig.tight_layout()
    out = os.path.join(FIGDIR, f"{name}.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


def main():
    sol_o, _ = lin.solve_ode()
    sol_d, _ = lin.solve_dae()
    _plot_pair("linear_qss", sol_o, sol_d, lin.biomass)

    sol_o, _ = quad.solve_ode()
    sol_d, _ = quad.solve_dae()
    _plot_pair("quadratic_qss", sol_o, sol_d, quad.biomass_positive_root)


if __name__ == "__main__":
    main()
