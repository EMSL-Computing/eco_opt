"""State/control, phase-portrait, and eco_opt vs analytical figures.

Run from the repository root::

    cd examples/microbial_growth
    python plot_microbial_growth.py

Output: figs/microbial_growth_states_control.png,
        figs/microbial_growth_phase_portrait.png,
        figs/microbial_growth_vs_analytical.png
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import microbial_growth as mg

HERE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(HERE, "figs")


def main():
    os.makedirs(FIG_DIR, exist_ok=True)

    sol, sys = mg.solve()
    assert sol.success, sol.return_status

    t = sol.time
    C = sol.x[:, 0]
    u = np.append(sol.u[:, 0], sol.u[-1, 0])
    G = mg.growth_rate(u)

    # -- Fig 1 equivalent: C(t), u(t), G(t) --------------------------------
    fig, axes = plt.subplots(3, 1, figsize=(6.5, 8), sharex=True)
    axes[0].plot(t, C, "o-", color="C0", ms=3)
    axes[0].set_ylabel("C(t) (resource pool)")
    axes[0].grid(alpha=0.3)
    axes[1].step(t, u, where="post", color="C1")
    axes[1].set_ylabel("u (uptake, control)")
    axes[1].grid(alpha=0.3)
    axes[2].plot(t, G, "-", color="C2")
    axes[2].set_ylabel("G (growth rate)")
    axes[2].set_xlabel("time")
    axes[2].grid(alpha=0.3)
    fig.suptitle(
        "Microbial growth OCP (Manzoni et al. 2022): eco_opt solution\n"
        f"J={sol.objective:.6f} (regression {mg.REGRESSION_OBJECTIVE:.6f}), "
        f"T={sol.tf:.6f}, C(T)={C[-1]:.2e}",
        fontsize=10,
    )
    fig.tight_layout()
    p1 = os.path.join(FIG_DIR, "microbial_growth_states_control.png")
    fig.savefig(p1, dpi=150)
    plt.close(fig)

    # -- Fig 2 equivalent: phase portraits (C-G, C-u, u-G) ------------------
    fig, axes = plt.subplots(3, 1, figsize=(5.5, 8))
    axes[0].plot(C, G, "-", color="C0")
    axes[0].set_xlabel("C"); axes[0].set_ylabel("G_opt"); axes[0].grid(alpha=0.3)
    axes[1].plot(C, u, "-", color="C1")
    axes[1].set_xlabel("C"); axes[1].set_ylabel("u_opt"); axes[1].grid(alpha=0.3)
    axes[2].plot(u, G, "-", color="C2")
    axes[2].set_xlabel("u_opt"); axes[2].set_ylabel("G_opt"); axes[2].grid(alpha=0.3)
    fig.suptitle("Microbial growth OCP: phase portraits (eco_opt)", fontsize=10)
    fig.tight_layout()
    p2 = os.path.join(FIG_DIR, "microbial_growth_phase_portrait.png")
    fig.savefig(p2, dpi=150)
    plt.close(fig)

    print("Wrote:", p1)
    print("Wrote:", p2)

    t_an, C_an, u_an, G_an, T_an, J_an = mg.analytical_solution()
    t_opt = sol.time
    C_opt = sol.x[:, 0]
    u_opt_t = sol.time_intervals
    u_opt = sol.u[:, 0]
    G_opt = mg.growth_rate(u_opt)
    J_an_t = np.concatenate(
        [[0.0], np.cumsum(0.5 * (G_an[1:] + G_an[:-1]) * np.diff(t_an))]
    )
    J_opt = np.concatenate([[0.0], np.cumsum(G_opt * np.diff(u_opt_t))])

    analytical_color = "#222222"
    state_color = "#0072B2"
    control_color = "#D55E00"
    secondary_color = "#009E73"
    fig, axes = plt.subplots(1, 3, figsize=(8.7, 3.0))

    ax = axes[0]
    ax.plot(t_an, u_an, color=analytical_color, lw=1.5, label="analytical")
    ax.step(
        u_opt_t[:-1],
        u_opt,
        where="post",
        color=control_color,
        lw=3,
        linestyle=":",
        label="eco_opt",
    )
    ax.set_xlabel("time")
    ax.set_ylabel("uptake $u$")
    ax.set_title("(a) Control")
    ax.legend(frameon=False, loc="upper right", handlelength=1.6)
    tf = float(sol.tf)
    ax.axvline(tf, color="0.4", lw=0.7, ls=":", zorder=0)
    ymin, ymax = ax.get_ylim()
    pad = 0.10 * (ymax - ymin)
    ax.set_ylim(ymin - pad, ymax)
    ax.annotate(
        r"$t_f$",
        xy=(tf, ymin),
        xytext=(tf - 0.5, ymin),
        ha="right",
        va="center",
        fontsize=12,
        color="0.25",
        arrowprops=dict(arrowstyle="<-", color="0.4", lw=0.7),
    )

    ax = axes[1]
    ax.plot(t_an, C_an, color=analytical_color, lw=1.5, label="analytical")
    ax.plot(t_opt, C_opt, color=state_color, lw=3, linestyle=":", label="eco_opt")
    ax.set_xlabel("time")
    ax.set_ylabel("substrate $C$")
    ax.set_title("(b) State")
    ax.legend(frameon=False, loc="upper right", handlelength=1.6)

    ax = axes[2]
    ax.plot(t_an, J_an_t, color=analytical_color, lw=1.5, label="analytical")
    ax.plot(
        u_opt_t,
        J_opt,
        color=secondary_color,
        lw=3,
        linestyle=":",
        label="eco_opt",
    )
    ax.set_xlabel("time")
    ax.set_ylabel("cumulative growth $J$")
    ax.set_title("(c) Objective")
    ax.legend(frameon=False, loc="upper left", handlelength=1.6)

    fig.tight_layout()
    p3 = os.path.join(FIG_DIR, "microbial_growth_vs_analytical.png")
    fig.savefig(p3, dpi=150)
    plt.close(fig)
    print("Wrote:", p3, f"(analytical T={T_an:.4f}, J={J_an:.4f})")


if __name__ == "__main__":
    main()
