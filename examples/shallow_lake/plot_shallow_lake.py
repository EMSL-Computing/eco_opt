"""Hysteresis of the lake, two management locals, and a restoration path."""
import os
import sys as _sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_sys.path.insert(0, os.path.dirname(__file__))
import shallow_lake as lake

FIGDIR = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIGDIR, exist_ok=True)


def _hysteresis_curve():
    P = np.linspace(1e-3, 2.6, 600)
    L = lake.loading_at_equilibrium(P)
    dL = np.gradient(L, P)
    folds = np.where(np.diff(np.sign(dL)))[0]
    return P, L, folds


def plot_all():
    P, L, folds = _hysteresis_curve()
    sol_lo, _ = lake.solve_management(L_guess=0.05)
    sol_hi, _ = lake.solve_management(L_guess=0.50)
    sol_r, _ = lake.solve_restoration()

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.6))

    ax = axes[0]
    ax.plot(L, P, color="0.25", lw=1.6)
    if folds.size >= 2:
        for i, lab in zip(folds[:2], ("fold (collapse)", "fold (restore)")):
            ax.plot(L[i], P[i], "o", color="C3", ms=5)
            ax.annotate(
                lab,
                xy=(L[i], P[i]),
                xytext=(8, 8 if i == folds[0] else -16),
                textcoords="offset points",
                fontsize=7,
                color="C3",
            )
    L_mark = 0.05
    eq = lake.equilibria_at_loading(L_mark)
    ax.axvline(L_mark, color="0.7", ls=":", lw=0.8)
    if eq.size >= 3:
        ax.plot(L_mark, eq[0], "o", color="C0", ms=6, label="clear")
        ax.plot(L_mark, eq[1], "s", color="0.5", ms=5, label="threshold")
        ax.plot(L_mark, eq[2], "o", color="C1", ms=6, label="turbid")
    ax.set_xlabel("constant loading $L$")
    ax.set_ylabel("phosphorus $P$")
    ax.set_title("(a) Ecological hysteresis")
    ax.legend(fontsize=7, loc="upper left")
    ax.set_xlim(0.0, 0.22)
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(sol_lo.time, sol_lo.x[:, 0], color="C1", label=r"$L_{\mathrm{guess}}\equiv 0.05$")
    ax.plot(sol_hi.time, sol_hi.x[:, 0], color="C0", label=r"$L_{\mathrm{guess}}\equiv 0.50$")
    ax.set_xlabel("time")
    ax.set_ylabel("$P(t)$")
    ax.set_title("(b) Two management locals")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    ax = axes[2]
    tL = sol_r.time_intervals[:-1]
    ax.step(tL, sol_r.u[:, 0], where="post", color="C2", label="$L(t)$")
    ax.axhline(lake.L_DESIRED, color="0.4", ls="--", lw=0.9, label=r"$L_{\mathrm{desired}}$")
    ax.set_xlabel("time")
    ax.set_ylabel("loading $L$")
    ax.set_title("(c) Restoration: cut, then ease")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    axb = ax.twinx()
    axb.plot(sol_r.time, sol_r.x[:, 0], color="C1", alpha=0.75)
    axb.set_ylabel("$P(t)$", color="C1")
    axb.axhline(lake.P_CLEAN, color="C1", ls=":", lw=0.8)

    fig.suptitle(
        "Shallow lake: recycling tipping point vs two optimal-control stories\n"
        f"management $J$={sol_lo.objective:.2f} (low guess) vs {sol_hi.objective:.2f} (high guess); "
        f"restore $P(T)$={sol_r.x[-1, 0]:.2f}",
        fontsize=10,
    )
    fig.tight_layout()
    out = os.path.join(FIGDIR, "shallow_lake.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    plot_all()
