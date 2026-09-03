"""Figures for the Cohen bang-bang and King–Roughgarden singular examples."""
import os
import sys as _sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_sys.path.insert(0, os.path.dirname(__file__))

import cohen_terminal_yield as cohen
import king_singular_allocation as king

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIGDIR, exist_ok=True)


def plot_cohen():
    sol, sys = cohen.solve()
    t_fine = np.linspace(0.0, cohen.T_SEASON, 400)
    veg_a, repro_a, u_a = cohen.analytic_trajectories(t_fine)
    tF = cohen.flowering_time_from_ocp(sol)

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.2), sharex=True)
    axes[0].plot(sol.time, sol.x[:, 0], "o", ms=3, color="C0", label="eco_opt veg $x_1$")
    axes[0].plot(sol.time, sol.x[:, 1], "s", ms=3, color="C1", label="eco_opt repro $x_2$")
    axes[0].plot(t_fine, veg_a, "-", color="C0", alpha=0.5, label="analytic $x_1$")
    axes[0].plot(t_fine, repro_a, "-", color="C1", alpha=0.5, label="analytic $x_2$")
    axes[0].axvline(cohen.ANALYTIC_SWITCH, color="k", ls="--", lw=1, label=f"analytic $t^*$={cohen.ANALYTIC_SWITCH:.2f}")
    axes[0].axvline(tF, color="0.5", ls=":", lw=1, label=f"OCP $t_F$={tF:.2f}")
    axes[0].set_ylabel("biomass")
    axes[0].legend(fontsize=8, loc="upper left")
    axes[0].grid(alpha=0.3)

    axes[1].step(sol.time_intervals[:-1], sol.u[:, 0], where="post", color="C2", label="eco_opt $u$")
    axes[1].plot(t_fine, u_a, "-", color="k", alpha=0.5, label="analytic $u$")
    axes[1].set_ylabel("allocation $u$ (vegetative fraction)")
    axes[1].set_xlabel("time")
    axes[1].set_ylim(-0.05, 1.05)
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)

    fig.suptitle(
        "Cohen (1971) terminal-yield allocation\n"
        f"J={sol.objective:.4f} (analytic {cohen.ANALYTIC_YIELD:.4f}), "
        f"$t^*$ analytic={cohen.ANALYTIC_SWITCH:.2f}",
        fontsize=11,
    )
    fig.tight_layout()
    out = os.path.join(FIGDIR, "cohen_terminal_yield.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


def plot_king():
    sol, sys = king.solve()
    t_fine = np.linspace(0.0, king.T_SEASON, 400)
    veg_a, repro_a, u_a = king.analytic_trajectories(t_fine)
    t2_ocp = king.t2_from_control(sol)
    # Ratio diagnostic is only meaningful on/near the singular arc; after t2
    # the formula x1/((T-1-t)x2) is not 1 and diverges as t → T-1.
    t_ratio_max = king.ANALYTIC_T2 + 0.35
    mask = sol.time <= t_ratio_max
    ratio = king.singular_ratio(sol.time[mask], sol.x[mask, 0], sol.x[mask, 1])
    t_ratio = t_fine[t_fine <= t_ratio_max]
    veg_r, repro_r, _ = king.analytic_trajectories(t_ratio)
    ratio_a = king.singular_ratio(t_ratio, veg_r, repro_r)

    fig, axes = plt.subplots(3, 1, figsize=(7.2, 9.0), sharex=True)
    axes[0].plot(sol.time, sol.x[:, 0], "o", ms=3, color="C0", label="eco_opt veg $x_1$")
    axes[0].plot(sol.time, sol.x[:, 1], "s", ms=3, color="C1", label="eco_opt repro $x_2$")
    axes[0].plot(t_fine, veg_a, "-", color="C0", alpha=0.5, label="analytic $x_1$")
    axes[0].plot(t_fine, repro_a, "-", color="C1", alpha=0.5, label="analytic $x_2$")
    axes[0].axvline(king.ANALYTIC_T2, color="k", ls="--", lw=1, label=f"analytic $t_2$={king.ANALYTIC_T2:.3f}")
    axes[0].set_ylabel("biomass")
    axes[0].legend(fontsize=8, loc="upper left")
    axes[0].grid(alpha=0.3)

    axes[1].plot(sol.time[mask], ratio, "o", ms=3, color="C3", label="eco_opt $x_1/((T-1-t)x_2)$")
    axes[1].plot(t_ratio, ratio_a, "-", color="k", alpha=0.5, label="analytic (=1 on singular arc)")
    axes[1].axhline(1.0, color="0.4", ls=":", lw=1)
    axes[1].axvline(king.ANALYTIC_T2, color="k", ls="--", lw=1)
    axes[1].set_ylabel("singular ratio")
    axes[1].set_ylim(-0.05, 1.6)
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)

    axes[2].step(sol.time_intervals[:-1], sol.u[:, 0], where="post", color="C2", label="eco_opt $u$")
    axes[2].plot(t_fine, u_a, "-", color="k", alpha=0.6, label="analytic $u^*$")
    axes[2].axvline(king.ANALYTIC_T2, color="k", ls="--", lw=1)
    axes[2].axvline(t2_ocp, color="0.5", ls=":", lw=1, label=f"OCP $t_2$={t2_ocp:.3f}")
    axes[2].set_ylabel("allocation $u$")
    axes[2].set_xlabel("time")
    axes[2].set_ylim(-0.05, 1.05)
    axes[2].legend(fontsize=8)
    axes[2].grid(alpha=0.3)

    term = sol.x[-1, 1] / sol.x[-1, 0]
    fig.suptitle(
        "King & Roughgarden (1982) graded allocation (start on singular manifold)\n"
        f"J={sol.objective:.4f} (analytic {king.ANALYTIC_OBJECTIVE:.4f}), "
        f"$x_2(T)/x_1(T)$={term:.3f} (analytic {king.ANALYTIC_TERMINAL_RATIO:.3f})",
        fontsize=11,
    )
    fig.tight_layout()
    out = os.path.join(FIGDIR, "king_singular_allocation.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    plot_cohen()
    plot_king()
