"""State/control figures for the isoperimetric-constraint example."""
import os
import sys as _sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_sys.path.insert(0, os.path.dirname(__file__))
import iso_constraint as iso

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIGDIR, exist_ok=True)


def _control_step(sol):
    t = np.asarray(sol.time_intervals)
    u = np.asarray(sol.u[:, 0])
    return t, np.append(u, u[-1])


def plot_iso_constraint():
    sol, _ = iso.solve()
    assert sol.success, sol.return_status
    err = abs(sol.objective - iso.REGRESSION_OBJECTIVE)
    note = (
        f"J={sol.objective:.6f} vs regression {iso.REGRESSION_OBJECTIVE:.6f} "
        f"(abs err {err:.2e}); x2(tf)={sol.x[-1, 1]:.6f} (target 10)"
    )

    fig, axes = plt.subplots(2, 1, figsize=(7, 5.8), sharex=True)
    axes[0].plot(sol.time, sol.x[:, 0], "o-", color="C0", ms=3, label="eco_opt")
    axes[0].set_ylabel("x1")
    axes[0].grid(alpha=0.3)
    axes[0].legend(loc="best", fontsize=8)
    axes[1].plot(sol.time, sol.x[:, 1], "o-", color="C0", ms=3)
    axes[1].set_ylabel("x2 (integral of u^2)")
    axes[1].set_xlabel("time")
    axes[1].grid(alpha=0.3)
    fig.suptitle(f"iso_constraint: states\n{note}", fontsize=10)
    fig.tight_layout()
    out = os.path.join(FIGDIR, "iso_constraint_states.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)

    fig, ax = plt.subplots(figsize=(7, 3.5))
    t, u = _control_step(sol)
    ax.step(t, u, where="post", color="C0", label="eco_opt")
    ax.set_xlabel("time")
    ax.set_ylabel("control u")
    ax.set_title("iso_constraint: control")
    ax.grid(alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    out = os.path.join(FIGDIR, "iso_constraint_control.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    plot_iso_constraint()
