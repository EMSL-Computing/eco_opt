"""Figures for the Iwasa–Roughgarden Liebig shoot/root/fruit example."""
import os
import sys as _sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_sys.path.insert(0, os.path.dirname(__file__))
import shoot_root as ir

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIGDIR, exist_ok=True)


def plot_shoot_root():
    sol, sys = ir.solve()
    t_fine = np.linspace(0.0, ir.T_SEASON, 500)
    S_a, R_a, F_a, uS_a, uR_a = ir.analytic_trajectories(t_fine)
    t_bal_ocp = ir.balance_time_from_ocp(sol)
    t_s_ocp = ir.fruit_switch_from_ocp(sol)

    fig, axes = plt.subplots(3, 1, figsize=(7.4, 9.2), sharex=True)

    axes[0].plot(sol.time, sol.x[:, 0], "o", ms=3, color="C2", label="eco_opt shoot $S$")
    axes[0].plot(sol.time, sol.x[:, 1], "s", ms=3, color="C3", label="eco_opt root $R$")
    axes[0].plot(sol.time, sol.x[:, 2], "^", ms=3, color="C1", label="eco_opt fruit $F$")
    axes[0].plot(t_fine, S_a, "-", color="C2", alpha=0.5, label="analytic $S$")
    axes[0].plot(t_fine, R_a, "-", color="C3", alpha=0.5, label="analytic $R$")
    axes[0].plot(t_fine, F_a, "-", color="C1", alpha=0.5, label="analytic $F$")
    axes[0].axvline(ir.ANALYTIC_T_BAL, color="0.2", ls="--", lw=1, label=f"$t_\\mathrm{{bal}}$={ir.ANALYTIC_T_BAL:.3f}")
    axes[0].axvline(ir.ANALYTIC_T_S, color="k", ls="--", lw=1, label=f"$t_s$={ir.ANALYTIC_T_S:.2f}")
    axes[0].set_ylabel("biomass")
    axes[0].legend(fontsize=7, loc="upper left", ncol=2)
    axes[0].grid(alpha=0.3)

    ratio = sol.x[:, 1] / np.maximum(sol.x[:, 0], 1e-12)
    ratio_a = R_a / np.maximum(S_a, 1e-12)
    axes[1].plot(sol.time, ratio, "o", ms=3, color="C4", label="eco_opt $R/S$")
    axes[1].plot(t_fine, ratio_a, "-", color="k", alpha=0.5, label="analytic $R/S$")
    axes[1].axhline(ir.ALPHA / ir.BETA, color="0.4", ls=":", label="balanced $α/β$")
    axes[1].axvline(ir.ANALYTIC_T_BAL, color="0.2", ls="--", lw=1)
    axes[1].axvline(ir.ANALYTIC_T_S, color="k", ls="--", lw=1)
    axes[1].set_ylabel("root/shoot ratio")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)

    t_u = sol.time_intervals[:-1]
    axes[2].step(t_u, sol.u[:, 0], where="post", color="C2", label="eco_opt $u_S$")
    axes[2].step(t_u, sol.u[:, 1], where="post", color="C3", label="eco_opt $u_R$")
    uF = 1.0 - sol.u[:, 0] - sol.u[:, 1]
    axes[2].step(t_u, uF, where="post", color="C1", label="eco_opt $u_F=1-u_S-u_R$")
    axes[2].plot(t_fine, uS_a, "-", color="C2", alpha=0.45)
    axes[2].plot(t_fine, uR_a, "-", color="C3", alpha=0.45)
    axes[2].plot(t_fine, 1.0 - uS_a - uR_a, "-", color="C1", alpha=0.45)
    axes[2].axvline(ir.ANALYTIC_T_BAL, color="0.2", ls="--", lw=1)
    axes[2].axvline(ir.ANALYTIC_T_S, color="k", ls="--", lw=1)
    axes[2].axvline(t_bal_ocp, color="0.5", ls=":", lw=1, label=f"OCP $t_\\mathrm{{bal}}$={t_bal_ocp:.2f}")
    axes[2].axvline(t_s_ocp, color="0.6", ls=":", lw=1, label=f"OCP $t_s$={t_s_ocp:.2f}")
    axes[2].set_ylabel("allocation")
    axes[2].set_xlabel("time")
    axes[2].set_ylim(-0.05, 1.05)
    axes[2].legend(fontsize=7, loc="center left", ncol=2)
    axes[2].grid(alpha=0.3)

    fig.suptitle(
        "Iwasa & Roughgarden (1984) three-phase shoot/root/fruit (Liebig $G=\\min(S,R)$)\n"
        f"F(T)={sol.x[-1, 2]:.4f} (analytic {ir.ANALYTIC_FRUIT:.4f}), "
        r"$t_{\mathrm{bal}}=\ln 5$, $t_s=T-1/\gamma=4$",
        fontsize=11,
    )
    fig.tight_layout()
    out = os.path.join(FIGDIR, "iwasa_roughgarden_shoot_root.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    plot_shoot_root()
