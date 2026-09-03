"""State/control figures for the nondimensional Goddard rocket example."""
import os
import sys as _sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_sys.path.insert(0, os.path.dirname(__file__))
import goddard_rocket_basic as grb

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIGDIR, exist_ok=True)


def _control_step(sol):
    t = np.asarray(sol.time_intervals)
    u = np.asarray(sol.u[:, 0])
    return t, np.append(u, u[-1])


def plot_goddard_rocket_basic():
    sol, _ = grb.solve()
    assert sol.success, sol.return_status
    err = abs(sol.objective - grb.REFERENCE_HEIGHT) / grb.REFERENCE_HEIGHT
    note = (
        f"h(T)={sol.objective:.6f} vs reference {grb.REFERENCE_HEIGHT} "
        f"(rel err {err:.2e}); T={sol.tf:.6f} vs reference {grb.REFERENCE_TF}"
    )
    labels = ["v (speed)", "h (height)", "m (mass)"]

    fig, axes = plt.subplots(3, 1, figsize=(7, 8.2), sharex=True)
    for i, ax in enumerate(axes):
        ax.plot(sol.time, sol.x[:, i], "o-", color="C0", ms=3, label="eco_opt" if i == 0 else None)
        ax.set_ylabel(labels[i])
        ax.grid(alpha=0.3)
    axes[0].legend(loc="best", fontsize=8)
    axes[-1].set_xlabel("time")
    fig.suptitle(f"goddard_rocket_basic: states\n{note}", fontsize=10)
    fig.tight_layout()
    out = os.path.join(FIGDIR, "goddard_rocket_basic_states.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)

    fig, ax = plt.subplots(figsize=(7, 3.5))
    t, u = _control_step(sol)
    ax.step(t, u, where="post", color="C0", label="eco_opt")
    ax.set_xlabel("time")
    ax.set_ylabel("T (thrust)")
    ax.set_title("goddard_rocket_basic: control")
    ax.grid(alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    out = os.path.join(FIGDIR, "goddard_rocket_basic_control.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    plot_goddard_rocket_basic()
