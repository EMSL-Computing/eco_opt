"""State/control and phase-portrait figures for the microbial growth example.

Run:  python examples/yop_microbial_growth/plot_microbial_growth.py
Output: figs/microbial_growth_states_control.png,
        figs/microbial_growth_phase_portrait.png
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


if __name__ == "__main__":
    main()
