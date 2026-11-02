"""Box 1 figure -- microbial uptake of a finite substrate.

(a) model schematic; (b) uptake, (c) substrate, (f) cumulative growth;
(d) shadow price of substrate; (e) uptake as a function of the shadow
price. Numerical optimal trajectories are compared with the analytical solution.

The optimum is the eco_opt solution, checked against the closed form of
Manzoni, Chakrawal and Ledder (2023). Requires eco_opt with costates
(``sol.costate``).

Run:
    python fig_box1_microbial.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import (  # noqa: E402
    apply_style, panel_label, save, shared_legend, TEXT_W,
    STATE_COLOR, CONTROL_COLOR, PRICE_COLOR, PAYOFF_COLOR, REF_COLOR,
    REF_LW, OPT_LW, OPT_LS,
)
from models import microbial as mg  # noqa: E402

apply_style()

A = mg.EMAX * mg.BETA
# Terminal uptake used in the analytical costate boundary condition.
U_TERMINAL = mg.R + np.sqrt(mg.R ** 2 + mg.R * mg.BETA)


def costate_closed_form(t, T):
    g_u = A * (mg.R + mg.BETA) / (U_TERMINAL + mg.BETA) ** 2
    lam_T = g_u / (1.0 - mg.MU * g_u)
    return lam_T * np.exp(mg.GAMMA * (np.asarray(t) - T))


def uptake_from_price(lam):
    return -mg.BETA + np.sqrt(A * (mg.R + mg.BETA) * (mg.MU + 1.0 / np.asarray(lam)))


def draw_schematic(ax):
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(0.08, 1.02)
    ax.axis("off")
    grey = "0.35"

    def pool(x, y, w, h, text, color):
        ax.add_patch(FancyBboxPatch((x, y), w, h, fc="white", ec=color, lw=1.0,
                                    boxstyle="round,pad=0.0,rounding_size=0.05"))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=7.5,
                linespacing=1.1)

    def arrow(p, q, color, rad=0.0, lw=0.9):
        ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=8,
                                     connectionstyle=f"arc3,rad={rad}", color=color, lw=lw))

    def label(x, y, text, color=grey, **kw):
        ax.text(x, y, text, fontsize=7, color=color, **kw)

    pool(-0.01, 0.38, 0.32, 0.30, "substrate\n$C$", STATE_COLOR)
    pool(0.69, 0.38, 0.32, 0.30, "microbes\n$B$", grey)
    # uptake (the investment) and recycling
    arrow((0.32, 0.57), (0.68, 0.57), CONTROL_COLOR, lw=1.5)
    label(0.50, 0.60, "uptake $u$", CONTROL_COLOR, ha="center", va="bottom")
    arrow((0.68, 0.48), (0.32, 0.48), grey)
    label(0.50, 0.45, r"recycling $\mu G$", ha="center", va="top")
    # respiration, loss, growth
    arrow((0.85, 0.69), (0.85, 0.86), grey)
    label(0.85, 0.88, r"CO$_2$", ha="center", va="bottom")
    arrow((0.15, 0.37), (0.15, 0.24), grey)
    label(0.15, 0.22, r"loss $\gamma C$", ha="center", va="top")
    arrow((0.85, 0.37), (0.85, 0.24), PAYOFF_COLOR)
    label(0.85, 0.22, "growth $G(u)$", PAYOFF_COLOR, ha="center", va="top")


def main():
    sol, _ = mg.solve()
    an = mg.analytical_solution()
    assert sol.success and sol.costate is not None

    T = float(sol.tf)
    t_u = sol.time_intervals
    u_opt = sol.u[:, 0]
    J_opt = np.concatenate([[0.0], np.cumsum(mg.growth_rate(u_opt) * np.diff(t_u))])
    t_lam, lam = sol.costate_time, sol.costate[:, 0]

    J_star = float(sol.objective)

    fig = plt.figure(figsize=(TEXT_W, 2.6))
    gs = fig.add_gridspec(2, 3, wspace=0.5, hspace=0.6, width_ratios=[1.6, 1, 1])
    ax = {k: fig.add_subplot(gs[i, j]) for k, (i, j) in
          zip("abcdef", [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)])}

    draw_schematic(ax["a"])
    panel_label(ax["a"], "a", x=0.02, y=0.98)

    # (b) uptake through time
    a_ = ax["b"]
    a_.plot(an["t"], an["u"], color=REF_COLOR, lw=REF_LW)
    a_.step(t_u[:-1], u_opt, where="post", color=CONTROL_COLOR, lw=OPT_LW, ls=OPT_LS)
    a_.set_ylim(0, 0.8)
    a_.set_xlim(0, 2.0)
    a_.set_xlabel("time")
    a_.set_ylabel("uptake $u$")
    panel_label(a_, "b")

    # (c) substrate
    a_ = ax["c"]
    a_.plot(an["t"], an["c"], color=REF_COLOR, lw=REF_LW)
    a_.plot(sol.time, sol.x[:, 0], color=STATE_COLOR, lw=OPT_LW, ls=OPT_LS, label="optimum")
    a_.set_xlim(0, 2.0)
    a_.set_ylim(0, 1.02)
    a_.set_xlabel("time")
    a_.set_ylabel(r"substrate $C/C_0$")
    panel_label(a_, "c")

    # (d) shadow price
    a_ = ax["d"]
    t_cf = np.linspace(0, T, 300)
    a_.plot(t_cf, costate_closed_form(t_cf, T), color=REF_COLOR, lw=REF_LW)
    a_.plot(t_lam, lam, color=PRICE_COLOR, lw=OPT_LW, ls=OPT_LS)
    a_.set_xlim(0, 2.0)
    a_.set_xlabel("time")
    a_.set_ylabel(r"shadow price $\lambda_C$")
    panel_label(a_, "d")

    # (e) uptake follows the price
    a_ = ax["e"]
    lam_grid = np.linspace(0.8 * lam.min(), 1.1 * lam.max(), 200)
    a_.plot(lam_grid, uptake_from_price(lam_grid), color=REF_COLOR, lw=REF_LW)
    a_.text(0.97, 0.97, "optimality\ncondition", transform=a_.transAxes,
            ha="right", va="top", fontsize=6.5)
    k = np.repeat(np.arange(len(u_opt)), sol.costate.shape[0] // len(u_opt))
    a_.plot(lam[::40], u_opt[k][::40], "o", ms=2.6, mfc=CONTROL_COLOR, mec="none")
    a_.annotate("$t=0$", xy=(lam[0], u_opt[0]), xytext=(lam[0] - 0.024, u_opt[0] - 0.19),
                fontsize=6, arrowprops=dict(arrowstyle="-", lw=0.5, color="0.4"))
    a_.annotate("$t=T$", xy=(lam[-1], u_opt[-1]), xytext=(lam[-1] - 0.01, u_opt[-1] + 0.14),
                fontsize=6, arrowprops=dict(arrowstyle="-", lw=0.5, color="0.4"))
    a_.set_xlabel(r"shadow price $\lambda_C$")
    a_.set_ylabel("uptake $u$")
    panel_label(a_, "e")

    # (f) cumulative growth
    a_ = ax["f"]
    J_an_t = np.concatenate(
        [[0.0], np.cumsum(0.5 * (an["g"][1:] + an["g"][:-1]) * np.diff(an["t"]))])
    a_.plot(an["t"], J_an_t, color=REF_COLOR, lw=REF_LW)
    a_.plot(t_u, J_opt, color=PAYOFF_COLOR, lw=OPT_LW, ls=OPT_LS)
    a_.set_xlim(0, 2.0)
    a_.set_ylim(0, 0.17)
    a_.set_xlabel("time")
    a_.set_ylabel("cumulative growth $J$")
    panel_label(a_, "f")

    dark = "0.2"
    shared_legend(fig, [
        ("analytical solution (Manzoni et al. 2023)", dict(color=REF_COLOR, lw=REF_LW)),
        ("numerical optimum", dict(color=dark, lw=OPT_LW, ls=OPT_LS)),
    ], y=0.95)
    save(fig, "fig_box1_microbial")

    # Numbers quoted in caption.md and results.md.
    g_T = lam[-1]
    print(f"T* = {T:.4f} (closed form {an['tf']:.4f});  J* = {J_star:.5f} (closed form {an['J']:.5f})")
    print(f"u(0) = {u_opt[0]:.3f}, u(T) = {u_opt[-1]:.3f}, terminal uptake = {U_TERMINAL:.3f}")
    print(f"lambda(0) = {lam[0]:.4f}, lambda(T) = {g_T:.4f}, ratio = {g_T / lam[0]:.3f}, "
          f"exp(gamma T) = {np.exp(mg.GAMMA * T):.3f}")
    err = np.max(np.abs(lam - costate_closed_form(t_lam, T))) / lam.max()
    print(f"max rel. error lambda = {err:.1e};  max |H| = {np.max(np.abs(sol.hamiltonian)):.1e}")



if __name__ == "__main__":
    main()
