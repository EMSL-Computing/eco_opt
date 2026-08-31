"""Reproduce (a Python/eco_opt analogue of) OptLignin's Figure 2: sensitivity
of the optimal decomposition trajectory to four parameters (vh_max, ro,
ARC0, mo), each swept over two values (base and 2x/3x), following
``examples/OptLignin/OptLignin1.0/Figure2.m``.

For each of 4 parameters x 2 values, solves the lignin decomposition OCP
(``lignin_decomposition.py``) and plots the resulting CT/CT(0), Co/Co(0)
(top row) and vo(t) (bottom row) trajectories, one column per swept
parameter -- mirroring the paper's 2x4 tiled-layout figure.

Run:  python examples/opt_lignin_figure2/plot_figure2.py
Output: figs/opt_lignin_figure2.png
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import lignin_decomposition as ld

HERE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(HERE, "figs")

# Sweep values, matching Figure2.m's ARC0/vh/ro/mo two-point sweeps.
VH = [ld.VH_LOOP, ld.VH_LOOP * 3]
RO = [ld.RO_LOOP, ld.RO_LOOP * 2]
ARC0 = [ld.ARC0_LOOP, ld.ARC0_LOOP * 2]
MO = [0.0, ld.MO_DEFAULT]

SWEEPS = [
    ("vh_max", VH, dict(vh_max=None)),
    ("ro", RO, dict(ro=None)),
    ("ARC0", ARC0, dict(ARC0=None)),
    ("mo", MO, dict(mo=None)),
]


def run_sweep(param_name, values):
    results = []
    for v in values:
        kwargs = dict(vh_max=ld.VH_LOOP, mo=ld.MO_DEFAULT, ro=ld.RO_LOOP, ARC0=ld.ARC0_LOOP)
        kwargs[param_name] = v
        sol, sys = ld.solve(**kwargs, control_intervals=30, polynomial_degree=3)
        results.append((v, sol))
    return results


def main():
    os.makedirs(FIG_DIR, exist_ok=True)

    fig, axes = plt.subplots(2, 4, figsize=(13, 5.5), sharex=True)
    linestyles = [":", "-"]
    titles = ["varying $v_H$", "varying $r_O$", "varying $ARC_0$", "varying $\\mu_O$"]

    for col, (param_name, values, _) in enumerate(SWEEPS):
        results = run_sweep(param_name, values)
        for (v, sol), ls in zip(results, linestyles):
            assert sol.success, f"{param_name}={v}: {sol.return_status}"
            t = sol.time
            CT = sol.x[:, 0]
            Co = sol.x[:, 1]
            vo = np.append(sol.u[:, 0], sol.u[-1, 0])
            axes[0, col].plot(t, CT / CT[0], ls, color="C0", lw=1.8)
            axes[0, col].plot(t, Co / Co[0], ls, color="C1", lw=1.8)
            axes[1, col].plot(t, vo, ls, color="k", lw=1.8, label=f"{param_name}={v:.4g}")
        axes[0, col].set_title(titles[col], fontsize=10)
        axes[0, col].grid(alpha=0.3)
        axes[1, col].grid(alpha=0.3)
        axes[1, col].set_xlabel("time [d]")
        axes[1, col].legend(fontsize=6, loc="best")

    axes[0, 0].set_ylabel("C remaining [-]")
    axes[1, 0].set_ylabel("$v_O$ [d$^{-1}$]")
    fig.suptitle(
        "OptLignin Figure 2 analogue (eco_opt): sensitivity of CT/Co/vo to "
        "vh_max, ro, ARC0, mo\n(solid = CT/CT0, dashed color = Co/Co0 top row; "
        "black = control vo bottom row)",
        fontsize=9,
    )
    fig.tight_layout()
    path = os.path.join(FIG_DIR, "opt_lignin_figure2.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("Wrote:", path)


if __name__ == "__main__":
    main()
