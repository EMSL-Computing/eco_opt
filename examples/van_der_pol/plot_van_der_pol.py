"""State/control figure for the Van der Pol example."""
import os
import sys as _sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

_sys.path.insert(0, os.path.dirname(__file__))
import van_der_pol as vdp

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIGDIR, exist_ok=True)


def plot_van_der_pol():
    sol, sys = vdp.solve()

    fig, axes = plt.subplots(3, 1, figsize=(7.0, 7.4), sharex=True)
    axes[0].plot(sol.time, sol.x[:, 0], "-", color="C0")
    axes[0].set_ylabel("$x_1$")
    axes[0].grid(alpha=0.3)
    axes[1].plot(sol.time, sol.x[:, 1], "-", color="C1")
    axes[1].set_ylabel("$x_2$")
    axes[1].grid(alpha=0.3)
    axes[2].step(sol.time_intervals[:-1], sol.u[:, 0], where="post", color="C2")
    axes[2].axhline(vdp.U_MIN, color="0.5", ls=":", lw=1)
    axes[2].axhline(vdp.U_MAX, color="0.5", ls=":", lw=1)
    axes[2].set_ylabel("$u$")
    axes[2].set_xlabel("time")
    axes[2].grid(alpha=0.3)

    fig.suptitle(
        "Van der Pol oscillator (CasADi collocation instance)\n"
        f"J={sol.objective:.6f} (regression {vdp.REGRESSION_OBJECTIVE:.6f}), "
        f"$x(T)$=({sol.x[-1, 0]:.3f}, {sol.x[-1, 1]:.3f})",
        fontsize=11,
    )
    fig.tight_layout()
    out = os.path.join(FIGDIR, "van_der_pol.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    plot_van_der_pol()
