"""State/control figures for the Bryson–Denham example."""
import os
import sys as _sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_sys.path.insert(0, os.path.dirname(__file__))
import bryson_denham as bd

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIGDIR, exist_ok=True)


def _control_step(sol):
    t = np.asarray(sol.time_intervals)
    u = np.asarray(sol.u[:, 0])
    return t, np.append(u, u[-1])


def plot_bryson_denham():
    sol, _ = bd.solve()
    assert sol.success, sol.return_status
    err = abs(sol.objective - bd.ANALYTIC_OBJECTIVE)
    note = (
        f"J={sol.objective:.6f} vs analytic {bd.ANALYTIC_OBJECTIVE} "
        f"(abs err {err:.2e}); max position={sol.x[:, 0].max():.6f} "
        f"(bound {bd.ANALYTIC_MAX_POSITION:.6f})"
    )

    fig, axes = plt.subplots(2, 1, figsize=(7, 5.8), sharex=True)
    axes[0].plot(sol.time, sol.x[:, 0], "o-", color="C0", ms=3, label="eco_opt")
    axes[0].axhline(bd.ANALYTIC_MAX_POSITION, color="k", ls="--", lw=1, label="x=1/9")
    axes[0].set_ylabel("position")
    axes[0].grid(alpha=0.3)
    axes[0].legend(loc="best", fontsize=8)
    axes[1].plot(sol.time, sol.x[:, 1], "o-", color="C0", ms=3)
    axes[1].set_ylabel("speed")
    axes[1].set_xlabel("time")
    axes[1].grid(alpha=0.3)
    fig.suptitle(f"bryson_denham: states\n{note}", fontsize=10)
    fig.tight_layout()
    out = os.path.join(FIGDIR, "bryson_denham_states.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)

    fig, ax = plt.subplots(figsize=(7, 3.5))
    t, u = _control_step(sol)
    ax.step(t, u, where="post", color="C0", label="eco_opt")
    ax.set_xlabel("time")
    ax.set_ylabel("u (acceleration)")
    ax.set_title("bryson_denham: control")
    ax.grid(alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    out = os.path.join(FIGDIR, "bryson_denham_control.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    plot_bryson_denham()
