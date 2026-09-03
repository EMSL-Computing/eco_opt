"""Figure for the two-phase Cohen example."""
import os
import sys as _sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_sys.path.insert(0, os.path.dirname(__file__))
import cohen_two_phase as cohen

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIGDIR, exist_ok=True)


def plot():
    sol, systems = cohen.solve()
    t_star = cohen.ANALYTIC_SWITCH
    t1 = np.linspace(0.0, t_star, 200)
    t2 = np.linspace(t_star, cohen.T_SEASON, 100)
    veg1 = cohen.X10 * np.exp(t1)
    repro1 = np.zeros_like(t1)
    veg2 = np.full_like(t2, cohen.X10 * np.exp(t_star))
    repro2 = veg2 * (t2 - t_star)

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.0), sharex=True)
    p1, p2 = sol.phases
    axes[0].plot(p1.time, p1.x[:, 0], "o", ms=3, color="C0", label="eco_opt veg (phase 1)")
    axes[0].plot(p2.time, p2.x[:, 0], "o", ms=3, color="C0", alpha=0.5, label="eco_opt veg (phase 2)")
    axes[0].plot(p1.time, p1.x[:, 1], "s", ms=3, color="C1", label="eco_opt repro (phase 1)")
    axes[0].plot(p2.time, p2.x[:, 1], "s", ms=3, color="C1", alpha=0.5, label="eco_opt repro (phase 2)")
    axes[0].plot(t1, veg1, "-", color="C0", alpha=0.45)
    axes[0].plot(t2, veg2, "-", color="C0", alpha=0.45)
    axes[0].plot(np.concatenate([t1, t2]), np.concatenate([repro1, repro2]), "-", color="C1", alpha=0.45)
    axes[0].axvline(t_star, color="k", ls="--", lw=1, label=f"analytic $t_s$={t_star:.2f}")
    axes[0].axvline(sol.switch_times[0], color="0.5", ls=":", lw=1)
    axes[0].set_ylabel("biomass")
    axes[0].legend(fontsize=8, loc="upper left")
    axes[0].grid(alpha=0.3)

    axes[1].axvspan(p1.t0, p1.tf, color="C2", alpha=0.15, label="phase 1: $\\dot x_1=x_1$")
    axes[1].axvspan(p2.t0, p2.tf, color="C3", alpha=0.15, label="phase 2: $\\dot x_2=x_1$")
    axes[1].axvline(sol.switch_times[0], color="k", ls="--", lw=1)
    axes[1].set_xlabel("time")
    axes[1].set_ylabel("phase")
    axes[1].set_yticks([])
    axes[1].legend(fontsize=8, loc="upper left")

    fig.suptitle(
        "Cohen (1971) two-phase transcription\n"
        f"$t_s$={sol.switch_times[0]:.5f} (analytic {t_star:.5f}), "
        f"$x_2(T)$={sol.objective:.4f} (analytic {cohen.ANALYTIC_YIELD:.4f})",
        fontsize=11,
    )
    fig.tight_layout()
    out = os.path.join(FIGDIR, "cohen_two_phase.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    plot()
