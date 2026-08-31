"""Reproduce Lindh et al. (2016) Fig. 2 (panels a-i): relative growth
rate, vegetative mass trajectories with optimal flowering time, and
optimal flowering time vs productivity, for exponential/logistic/WBE
growth. See README.md and optimal_flowering.py for the underlying model.

Uses the fast 1-D (closed-form ODE + scalar optimization) approach from
optimal_flowering.py to sweep productivity densely; the OCP-based
bang-bang verification (agreeing with the 1-D result at a few sample
points) lives in optimal_flowering.py / test_lindh_2016_flowering.py.
"""
import os
import sys as _sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_sys.path.insert(0, os.path.dirname(__file__))
from optimal_flowering import (
    GROWTH_MODELS, V0, VMAX, SEASON_END,
    integrate_V, optimal_flowering_time_1d,
)

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIGDIR, exist_ok=True)

GROWTH_ORDER = ["exponential", "logistic", "wbe"]
GROWTH_TITLES = {"exponential": "Exponential growth", "logistic": "Logistic growth", "wbe": "West-Brown-Enquist growth"}


def relative_growth_rate(growth, V, P=1.0, vmax=VMAX):
    U_fn = GROWTH_MODELS[growth]
    Vc = np.maximum(V, 1e-9)
    F = P * U_fn(Vc, vmax) if growth != "exponential" else P * U_fn(Vc)
    return F / Vc


def vegetative_trajectory(growth, P, tF, n=200, v0=V0, vmax=VMAX):
    """Sample V(t) for t in [0, tF] (uncontrolled ODE, valid pre-flowering)."""
    from scipy.integrate import solve_ivp
    U_fn = GROWTH_MODELS[growth]

    def rhs(t, V):
        Vc = max(V[0], 1e-9)
        return [P * U_fn(Vc, vmax) if growth != "exponential" else P * U_fn(Vc)]

    t_eval = np.linspace(0, tF, n)
    sol = solve_ivp(rhs, [0, tF], [v0], t_eval=t_eval, rtol=1e-9, atol=1e-12)
    return sol.t, sol.y[0]


def make_figure2(P_levels=(0.5, 1.0, 2.0, 4.0, 8.0), P_sweep=None):
    if P_sweep is None:
        P_sweep = np.linspace(0.3, 15.0, 80)

    fig, axes = plt.subplots(3, 3, figsize=(11, 9))
    cmap = matplotlib.colormaps["viridis"]
    colors = [cmap(x) for x in np.linspace(0, 1, len(P_levels))]

    # Row 1 (a-c): relative growth rate vs V
    V_grid = np.linspace(1e-3, VMAX * 0.999, 200)
    for col, growth in enumerate(GROWTH_ORDER):
        ax = axes[0, col]
        for i, P in enumerate(P_levels):
            rgr = relative_growth_rate(growth, V_grid, P=1.0)  # RGR indep. of P for these models (F=P*U(V), F/V=P*U(V)/V)
            # RGR scales linearly with P; show P-independent SHAPE curves at
            # representative P (paper shows a family of curves for different P
            # in (b),(c) since RGR = P * U(V)/V does depend on P).
            ax.plot(V_grid, P * (relative_growth_rate(growth, V_grid, P=1.0)), color=colors[i], lw=1.5)
        ax.set_title(f"({chr(97+col)})", loc="left")
        ax.set_xlabel("Vegetative mass (V)")
        if col == 0:
            ax.set_ylabel("Relative growth rate")
        ax.set_xlim(0, VMAX)

    # Row 2 (d-f): vegetative mass trajectories + optimal flowering time markers
    for col, growth in enumerate(GROWTH_ORDER):
        ax = axes[1, col]
        for i, P in enumerate(P_levels):
            tF = optimal_flowering_time_1d(growth, P)
            t, V = vegetative_trajectory(growth, P, SEASON_END)  # full-season "had flowering not occurred"
            ax.plot(t, V, color=colors[i], lw=1, ls=":", alpha=0.6)
            t_f_traj, V_f_traj = vegetative_trajectory(growth, P, tF)
            ax.plot(t_f_traj, V_f_traj, color=colors[i], lw=1.8)
            ax.plot(tF, V_f_traj[-1], "o", color=colors[i], ms=5)
        ax.set_title(f"({chr(100+col)})", loc="left")
        ax.set_xlabel("Time (t)")
        if col == 0:
            ax.set_ylabel("Vegetative mass (V)")
        ax.set_xlim(0, SEASON_END)

    # Row 3 (g-i): optimal flowering time vs productivity -- THE key result
    for col, growth in enumerate(GROWTH_ORDER):
        ax = axes[2, col]
        tF_vals = np.array([optimal_flowering_time_1d(growth, P) for P in P_sweep])
        ax.plot(P_sweep, tF_vals, color="black", lw=2)
        for i, P in enumerate(P_levels):
            tF = optimal_flowering_time_1d(growth, P)
            ax.plot(P, tF, "o", color=colors[i], ms=6)
        ax.set_title(f"({chr(103+col)})", loc="left")
        ax.set_xlabel("Productivity (P)")
        if col == 0:
            ax.set_ylabel("Optimal flowering time (t*_F)")
        ax.set_ylim(0, SEASON_END)

    for ax in axes.flat:
        ax.grid(True, alpha=0.4)

    fig.suptitle("Lindh et al. (2016) Fig. 2 reproduction: exponential vs logistic vs WBE growth\n"
                 "(row 3: constrained growth flips the sign of d(optimal flowering time)/dP)")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(FIGDIR, "figure2_reproduction.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)
    return P_sweep, {g: np.array([optimal_flowering_time_1d(g, P) for P in P_sweep]) for g in GROWTH_ORDER}


if __name__ == "__main__":
    P_sweep, tF_curves = make_figure2()
    for g in GROWTH_ORDER:
        tf = tF_curves[g]
        is_monotone_increasing = np.all(np.diff(tf) >= -1e-9)
        argmax_idx = np.argmax(tf)
        interior_max = 0 < argmax_idx < len(tf) - 1
        print(f"{g:12s}: monotone increasing={is_monotone_increasing} "
              f"interior max at P={P_sweep[argmax_idx]:.3f} (idx {argmax_idx}/{len(tf)-1})")
