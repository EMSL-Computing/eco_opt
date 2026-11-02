"""Box 2 figure -- ligninolytic investment during litter decomposition.

(a) model schematic; (b) optimal oxidative investment v_O(t) at two costs
r_O; (c) remaining non-aromatic and aromatic
carbon; (d) shadow prices of non-aromatic and aromatic carbon; (e)
marginal return to initiating ligninolysis, dH/dv_O at v_O = 0, whose sign
change marks the onset of investment; (f) onset (lag) time across r_O.

Model and parameters from Chakrawal, Lindahl and Manzoni (2024).
Requires eco_opt with costates (``sol.costate``).

Run:
    python fig_box2_lignin.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.path import Path as MplPath

sys.path.insert(0, str(Path(__file__).resolve().parent))
from style import (  # noqa: E402
    apply_style, panel_label, save, shared_legend, TEXT_W,
    STATE_COLOR, CONTROL_COLOR, PRICE_COLOR, PAYOFF_COLOR,
)
from models import lignin as ld  # noqa: E402

apply_style()

RO_LOW, RO_HIGH = 25.0, 50.0
RO_SWEEP = np.array([5.0, 8.0, 12.0, 16.0, 20.0, 25.0, 35.0, 50.0, 75.0, 100.0])
N_INTERVALS = 100
IPOPT = {"max_iter": 3000, "tol": 1e-10}

STATE_LIGHT = "#56B4E9"
PRICE_DARK = "#882255"
LS_LOW = "-"
LS_HIGH = (0, (1.2, 1.2))
LW = 1.5


def g_kernel(L):
    return np.exp(-(L / ld.A) ** ld.B)


def lag_time(t_edges, vo, frac=0.05):
    """First time v_O exceeds a fraction of its maximum (source definition)."""
    hit = np.where(vo > frac * np.max(vo))[0]
    return float(t_edges[hit[0]]) if hit.size else np.nan


def solve(ocp, sys_, ro):
    ld.apply_scenario(ocp, sys_, vh_max=ld.VH_LOOP, mo=ld.MO_DEFAULT, ro=ro)
    sol = ocp.optimize(ipopt_options=IPOPT)
    assert sol.success, f"r_O={ro}: {sol.return_status}"
    t = sol.costate_time
    lam = sol.costate
    CT = np.interp(t, sol.time, sol.x[:, 0])
    Co = np.interp(t, sol.time, sol.x[:, 1])
    Dh = ld.VH_LOOP * g_kernel(Co / CT) * (CT - Co)
    # With states (C_T, C_O): non-aromatic C_H = C_T - C_O, so
    # dJ/dC_H = lambda_T and dJ/dC_O (at fixed C_H) = lambda_T + lambda_O.
    lam_H = lam[:, 0]
    lam_A = lam[:, 0] + lam[:, 1]
    # dH/dv_O at v_O = 0: growth gained from oxidation, valued at
    # 1 + lambda_T + m_O lambda_O, minus the efficiency cost on all
    # decomposition, minus the value of the aromatic carbon consumed.
    w = 1.0 + lam[:, 0] + ld.MO_DEFAULT * lam[:, 1]
    dH0 = w * (ld.EMAX * Co - ro * Dh) - lam_A * Co
    return {
        "sol": sol, "J": float(sol.objective), "t": sol.time, "t_u": sol.time_intervals,
        "vo": sol.u[:, 0], "CH": sol.x[:, 0] - sol.x[:, 1], "CA": sol.x[:, 1],
        "t_lam": t, "lam_H": lam_H, "lam_A": lam_A, "dH0": dH0,
        "tau": lag_time(sol.time_intervals, sol.u[:, 0]),
    }


def draw_schematic(ax):
    ax.set_xlim(-0.02, 1.30)
    ax.set_ylim(-0.38, 1.38)
    ax.axis("off")
    grey = "0.35"

    def pool(x, y, w, h, text, color):
        ax.add_patch(FancyBboxPatch((x, y), w, h, fc="white", ec=color, lw=1.0,
                                    boxstyle="round,pad=0.0,rounding_size=0.035"))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=7,
                linespacing=1.2)

    def arrow(p, q, color, rad=0.0, lw=0.9):
        ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=8,
                                     connectionstyle=f"arc3,rad={rad}", color=color, lw=lw))

    def curve(p, c, q, color, lw=0.9):
        # Quadratic Bezier in data coordinates, so the bulge does not depend on the axes aspect.
        path = MplPath([p, c, q], [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3])
        ax.add_patch(FancyArrowPatch(path=path, arrowstyle="-|>", mutation_scale=8,
                                     color=color, lw=lw))

    def label(x, y, text, color=grey, **kw):
        ax.text(x, y, text, fontsize=6.5, color=color, **kw)

    pool(0.00, 0.72, 0.50, 0.36, "non-aromatic\n$C_H$", STATE_COLOR)
    pool(0.00, 0.00, 0.50, 0.36, "aromatic\n$C_O$", STATE_LIGHT)
    pool(0.80, 0.39, 0.34, 0.27, "microbes", grey)
    # Uptake: hydrolysis labelled above its arrow, ligninolysis below its arrow.
    arrow((0.51, 0.86), (0.79, 0.61), grey)
    label(0.56, 0.86, "hydrolysis", ha="left", va="bottom")
    arrow((0.51, 0.20), (0.79, 0.45), CONTROL_COLOR, lw=1.5)
    label(0.53, 0.17, "ligninolysis\n" + r"$v_O C_O$", CONTROL_COLOR,
          ha="left", va="top", linespacing=1.2)
    # Recycling routed around the outside of the pools.
    curve((1.10, 0.67), (1.15, 1.35), (0.20, 1.10), grey)
    label(0.70, 1.22, r"$(1-\mu_O)G$", ha="center", va="bottom")
    curve((1.10, 0.38), (1.15, -0.42), (0.20, -0.01), grey)
    label(0.70, -0.21, r"necromass $\mu_O G$", ha="center", va="top")
    arrow((1.15, 0.525), (1.29, 0.525), grey)
    label(1.22, 0.60, r"CO$_2$", ha="center", va="bottom")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot-data", type=Path,
                        help="Save/reuse numerical trajectories during figure layout edits.")
    args = parser.parse_args()
    keys = ("J", "t", "t_u", "vo", "CH", "CA", "t_lam", "lam_H", "lam_A", "dH0", "tau")
    if args.plot_data is not None and args.plot_data.exists():
        with np.load(args.plot_data) as data:
            low = {k: data[f"low_{k}"] for k in keys}
            high = {k: data[f"high_{k}"] for k in keys}
            taus = data["taus"]
    else:
        ocp, sys_ = ld.compile_ocp(control_intervals=N_INTERVALS)
        low = solve(ocp, sys_, RO_LOW)
        high = solve(ocp, sys_, RO_HIGH)
        taus = np.array([solve(ocp, sys_, ro)["tau"] for ro in RO_SWEEP])
        if args.plot_data is not None:
            records = {f"{name}_{k}": rec[k] for name, rec in (("low", low), ("high", high))
                       for k in keys}
            np.savez(args.plot_data, taus=taus, **records)

    fig = plt.figure(figsize=(TEXT_W, 4.05))
    gs = fig.add_gridspec(2, 3, wspace=0.5, hspace=0.6, width_ratios=[2.0, 1, 1],
                          height_ratios=[1.45, 1])
    ax = {k: fig.add_subplot(gs[i, j]) for k, (i, j) in
          zip("abcdef", [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)])}

    draw_schematic(ax["a"])
    panel_label(ax["a"], "a", x=0.02, y=0.98)

    # (b) investment through time
    a_ = ax["b"]
    for rec, ls in ((low, LS_LOW), (high, LS_HIGH)):
        a_.step(rec["t_u"][:-1], 1e3 * rec["vo"], where="post", color=CONTROL_COLOR,
                ls=ls, lw=LW)
    a_.annotate(r"$\tau$", xy=(low["tau"], 0), xytext=(low["tau"] - 260, 0.35),
                fontsize=7, arrowprops=dict(arrowstyle="-", lw=0.5, color="0.4"))
    a_.set_xlim(0, ld.TT_DEFAULT)
    a_.set_ylim(0, 2.3)
    a_.set_xlabel("time (d)")
    a_.set_ylabel(r"ligninolytic rate $v_O$ ($10^{-3}$ d$^{-1}$)")
    panel_label(a_, "b")

    # (c) remaining carbon (r_O = 25)
    a_ = ax["c"]
    a_.plot(low["t"], low["CH"] / low["CH"][0], color=STATE_COLOR, lw=LW)
    a_.plot(low["t"], low["CA"] / low["CA"][0], color=STATE_LIGHT, lw=LW)
    a_.text(60, 0.10, "non-aromatic", color=STATE_COLOR, fontsize=6, ha="left")
    a_.text(420, 1.24, "aromatic", color=STATE_LIGHT, fontsize=6, ha="left")
    a_.set_xlim(0, ld.TT_DEFAULT)
    a_.set_ylim(0, 1.35)
    a_.set_xlabel("time (d)")
    a_.set_ylabel("remaining C (relative)")
    panel_label(a_, "c")

    # (d) shadow prices (r_O = 25)
    a_ = ax["d"]
    a_.axhline(0, color="0.6", lw=0.6)
    a_.plot(low["t_lam"], low["lam_H"], color=PRICE_COLOR, lw=LW)
    a_.plot(low["t_lam"], low["lam_A"], color=PRICE_DARK, lw=LW)
    a_.text(60, 0.37, r"non-aromatic $\lambda_H$", color=PRICE_COLOR, fontsize=6, va="bottom")
    a_.text(1150, -0.34, r"aromatic $\lambda_O$", color=PRICE_DARK, fontsize=6, va="top")
    a_.set_xlim(0, ld.TT_DEFAULT)
    a_.set_ylim(-0.42, 0.48)
    a_.set_yticks([-0.25, 0, 0.25])
    a_.set_xlabel("time (d)")
    a_.set_ylabel("shadow price (per unit C)")
    panel_label(a_, "d")

    # (e) marginal return to initiating ligninolysis
    a_ = ax["e"]
    a_.axhline(0, color="0.6", lw=0.6)
    for rec, ls in ((low, LS_LOW), (high, LS_HIGH)):
        a_.plot(rec["t_lam"], 1e3 * rec["dH0"], color=PAYOFF_COLOR, ls=ls, lw=LW)
    a_.text(0.03, 0.97, r"$>0$: invest", transform=a_.transAxes, ha="left", va="top",
            fontsize=6, color="0.4")
    a_.text(0.97, 0.04, r"$<0$: wait", transform=a_.transAxes, ha="right", va="bottom",
            fontsize=6, color="0.4")
    a_.set_xlim(0, ld.TT_DEFAULT)
    a_.set_ylim(None, 7)
    a_.set_xlabel("time (d)")
    a_.set_ylabel("marginal return to\ninitiation ($10^{-3}$)")
    panel_label(a_, "e")

    # (f) onset of investment across costs
    a_ = ax["f"]
    a_.plot(RO_SWEEP, taus, "o-", color=CONTROL_COLOR, ms=3, lw=1.0)
    for ro, rec in ((RO_LOW, low), (RO_HIGH, high)):
        a_.plot(ro, rec["tau"], "o", ms=5, mfc="none", mec="0.2", mew=0.8)
    a_.set_xscale("log")
    a_.set_xlabel(r"cost of ligninolysis $r_O$ (d)")
    a_.set_ylabel(r"lag of ligninolysis $\tau$ (d)")
    panel_label(a_, "f")

    shared_legend(fig, [
        (rf"optimum, $r_O={RO_LOW:.0f}$ d", dict(color="0.2", lw=LW, ls=LS_LOW)),
        (rf"optimum, $r_O={RO_HIGH:.0f}$ d", dict(color="0.2", lw=LW, ls=LS_HIGH)),
    ], y=0.95)
    save(fig, "fig_box2_lignin")

    # Numbers quoted in caption.md and results.md.
    for name, rec in (("r_O=25", low), ("r_O=50", high)):
        s = np.sign(rec["dH0"])
        cross = rec["t_lam"][np.where(np.diff(s) != 0)[0]]
        print(f"{name}: J* = {rec['J']:.5f}, tau (5% rule) = {rec['tau']:.0f} d, "
              f"dH/dv_O|0 changes sign at {np.round(cross).tolist()} d, "
              f"max v_O = {rec['vo'].max():.2e}")
        print(f"   lambda_H(0) = {rec['lam_H'][0]:.3f}, lambda_O(0) = {rec['lam_A'][0]:.3f}, "
              f"min lambda_O = {rec['lam_A'].min():.3f}")
        print(f"   C_H(T)/C_H0 = {rec['CH'][-1] / rec['CH'][0]:.3f}, "
              f"C_O(T)/C_O0 = {rec['CA'][-1] / rec['CA'][0]:.3f}")
    print("tau(r_O):", dict(zip(RO_SWEEP.tolist(), np.round(taus).tolist())))


if __name__ == "__main__":
    main()
