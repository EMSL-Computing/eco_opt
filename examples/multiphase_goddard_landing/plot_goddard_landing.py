"""State, control, and drag-force figures for Goddard landing."""
import os
import sys as _sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

_sys.path.insert(0, os.path.dirname(__file__))
import goddard_landing as gl

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIGDIR, exist_ok=True)


def plot():
    sol, systems = gl.solve()
    p1, p2 = sol.phases
    sys1, sys2 = systems

    fig, axes = plt.subplots(4, 1, figsize=(7.2, 9.0), sharex=True)
    labels = [("velocity (m/s)", 0), ("height (m)", 1), ("mass (kg)", 2)]
    for ax, (ylab, i) in zip(axes[:3], labels):
        ax.plot(p1.time, p1.x[:, i], "-", color="C0", label="ascent")
        ax.plot(p2.time, p2.x[:, i], "-", color="C1", label="landing")
        ax.axvline(sol.switch_times[0], color="k", ls="--", lw=1)
        ax.set_ylabel(ylab)
        ax.grid(alpha=0.3)
        if i == 0:
            ax.legend(fontsize=8)
    axes[3].step(p1.time_intervals[:-1], p1.u[:, 0], where="post", color="C0")
    axes[3].step(p2.time_intervals[:-1], p2.u[:, 0], where="post", color="C1")
    axes[3].axvline(sol.switch_times[0], color="k", ls="--", lw=1)
    axes[3].set_ylabel("F (kg/s)")
    axes[3].set_xlabel("time (s)")
    axes[3].grid(alpha=0.3)
    fig.suptitle(
        "Goddard landing (two-phase)\n"
        f"peak $h$={sol.objective:.0f} m at $t_s$={sol.switch_times[0]:.1f} s, "
        f"$t_f$={sol.tf:.1f} s, $h(T)$={sol.x[-1, 1]:.2e} m",
        fontsize=11,
    )
    fig.tight_layout()
    out1 = os.path.join(FIGDIR, "goddard_landing_states.png")
    fig.savefig(out1, dpi=150)
    plt.close(fig)

    drag1 = p1.signal(sys1.y.drag.force)
    drag2 = p2.signal(sys2.y.drag.force)
    fig = plt.figure(figsize=(6.5, 5.5))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(p1.x[:, 1], p1.x[:, 0], drag1, color="C0", label="ascent")
    ax.plot(p2.x[:, 1], p2.x[:, 0], drag2, color="C1", label="landing")
    ax.set_xlabel("h (m)")
    ax.set_ylabel("v (m/s)")
    ax.set_zlabel("drag (N)")
    ax.legend(fontsize=8)
    fig.suptitle("Goddard landing: drag vs height and velocity", fontsize=11)
    fig.tight_layout()
    out2 = os.path.join(FIGDIR, "goddard_landing_drag3d.png")
    fig.savefig(out2, dpi=150)
    plt.close(fig)
    print("wrote", out1)
    print("wrote", out2)


if __name__ == "__main__":
    plot()
