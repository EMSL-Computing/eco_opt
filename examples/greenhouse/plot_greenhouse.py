"""States, heating, and exogenous light / outside temperature."""
import os
import sys as _sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_sys.path.insert(0, os.path.dirname(__file__))
import greenhouse as gh

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIGDIR, exist_ok=True)


def plot_greenhouse():
    sol, sys = gh.solve()
    t = sol.time
    t_u = sol.time_intervals[:-1]
    I = np.array([float(gh.light(ti)) for ti in t])
    T_out = np.array([float(gh.outside_temp(ti)) for ti in t])

    fig, axes = plt.subplots(3, 1, figsize=(7.2, 7.6), sharex=True)
    axes[0].plot(t, sol.x[:, 0] * 1200.0, "-", color="C0")
    axes[0].set_ylabel(r"$1200\,x_1$ (crop)")
    axes[0].grid(alpha=0.3)
    axes[1].plot(t, sol.x[:, 1], "-", color="C1", label="indoor")
    axes[1].plot(t, T_out, "--", color="0.4", label="outside")
    axes[1].set_ylabel(r"temperature [$^\circ$C]")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)
    axes[2].plot(t, I / 40.0, "-", color="C2", label=r"$I/40$")
    axes[2].step(t_u, sol.u[:, 0], where="post", color="C3", label="heat $u$")
    axes[2].set_ylabel(r"$I/40$, $u$")
    axes[2].set_xlabel("time [h]")
    axes[2].legend(fontsize=8)
    axes[2].grid(alpha=0.3)
    fig.suptitle(
        "Greenhouse climate control (interpolated light and $T_{out}$)\n"
        f"J={sol.objective:.4f} (regression {gh.REGRESSION_OBJECTIVE:.4f}), "
        f"crop$(T)$={sol.x[-1, 0]:.4f} kg m$^{{-2}}$",
        fontsize=11,
    )
    fig.tight_layout()
    out = os.path.join(FIGDIR, "greenhouse.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    plot_greenhouse()
