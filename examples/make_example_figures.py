"""Generate one state/control comparison figure for each validated eco_opt
example (bang_free_tf, bryson_denham, goddard_rocket_basic, iso_constraint,
non_diff), on top of the existing dedicated Goddard free-tf plot script.

Each figure overlays the eco_opt trajectory against whatever reference is
available for that example (analytic value / literature reference /
regression-tracked own value -- see each example module's docstring), and
annotates the objective error so the plot is self-documenting.

This script does not modify any existing example or package file; it only
imports the existing examples/*.py modules and plots their solve() output.

Run:  python examples/make_example_figures.py
Output: examples/figs/<name>_states.png, examples/figs/<name>_control.png
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import bang_free_tf as bft
import bryson_denham as bd
import goddard_rocket_basic as grb
import iso_constraint as iso
import non_diff as nd

HERE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(HERE, "figs")


def _control_step(sol):
    """Return (t, u) step arrays for a piecewise-constant control."""
    t = np.asarray(sol.time_intervals)
    u = np.asarray(sol.u[:, 0])
    return t, np.append(u, u[-1])


def plot_example(name, sol, state_labels, ref_note, control_label="control u"):
    os.makedirs(FIG_DIR, exist_ok=True)
    n_states = sol.x.shape[1]

    fig, axes = plt.subplots(n_states, 1, figsize=(7, 2.4 * n_states + 1), sharex=True)
    if n_states == 1:
        axes = [axes]
    for i, ax in enumerate(axes):
        ax.plot(sol.time, sol.x[:, i], "o-", color="C0", ms=3, label="eco_opt")
        ax.set_ylabel(state_labels[i])
        ax.grid(alpha=0.3)
    axes[0].legend(loc="best", fontsize=8)
    axes[-1].set_xlabel("time")
    fig.suptitle(f"{name}: states\n{ref_note}", fontsize=10)
    fig.tight_layout()
    states_path = os.path.join(FIG_DIR, f"{name}_states.png")
    fig.savefig(states_path, dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 3.5))
    t, u = _control_step(sol)
    ax.step(t, u, where="post", color="C0", label="eco_opt")
    ax.set_xlabel("time")
    ax.set_ylabel(control_label)
    ax.set_title(f"{name}: control")
    ax.grid(alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    control_path = os.path.join(FIG_DIR, f"{name}_control.png")
    fig.savefig(control_path, dpi=150)
    plt.close(fig)

    return states_path, control_path


def main():
    written = []

    # -- bang_free_tf: analytic bang-bang double integrator ---------------
    sol, _ = bft.solve()
    assert sol.success, sol.return_status
    err = abs(sol.tf - bft.ANALYTIC_TF)
    note = f"tf={sol.tf:.6f} vs analytic {bft.ANALYTIC_TF} (abs err {err:.2e})"
    written += plot_example("bang_free_tf", sol, ["x1 (position)", "x2 (velocity)"],
                             note, control_label="u (accel)")

    # -- bryson_denham: analytic closed-form objective ---------------------
    sol, _ = bd.solve()
    assert sol.success, sol.return_status
    err = abs(sol.objective - bd.ANALYTIC_OBJECTIVE)
    note = (f"J={sol.objective:.6f} vs analytic {bd.ANALYTIC_OBJECTIVE} (abs err {err:.2e}); "
            f"max position={sol.x[:,0].max():.6f} (bound {bd.ANALYTIC_MAX_POSITION:.6f})")
    written += plot_example("bryson_denham", sol, ["position", "speed"], note,
                             control_label="u (acceleration)")

    # -- goddard_rocket_basic: literature (Petit et al. 2008) reference ---
    sol, _ = grb.solve()
    assert sol.success, sol.return_status
    err = abs(sol.objective - grb.REFERENCE_HEIGHT) / grb.REFERENCE_HEIGHT
    note = (f"h(T)={sol.objective:.6f} vs reference {grb.REFERENCE_HEIGHT} (rel err {err:.2e}); "
            f"T={sol.tf:.6f} vs reference {grb.REFERENCE_TF}")
    written += plot_example("goddard_rocket_basic", sol,
                             ["v (speed)", "h (height)", "m (mass)"], note,
                             control_label="T (thrust)")

    # -- iso_constraint: regression-tracked own value ----------------------
    sol, _ = iso.solve()
    assert sol.success, sol.return_status
    err = abs(sol.objective - iso.REGRESSION_OBJECTIVE)
    note = (f"J={sol.objective:.6f} vs regression value {iso.REGRESSION_OBJECTIVE:.6f} "
            f"(abs err {err:.2e}); x2(tf)={sol.x[-1,1]:.6f} (isoperimetric target 10)")
    written += plot_example("iso_constraint", sol, ["x1", "x2 (integral of u^2)"], note)

    # -- non_diff: regression-tracked own value -----------------------------
    sol, _ = nd.solve()
    assert sol.success, sol.return_status
    err = abs(sol.objective - nd.REGRESSION_OBJECTIVE)
    note = (f"J=x3(tf)={sol.objective:.6f} vs regression value "
            f"{nd.REGRESSION_OBJECTIVE:.6f} (abs err {err:.2e})")
    written += plot_example("non_diff", sol, ["x1", "x2", "x3 (running cost)"], note)

    for p in written:
        print("Wrote:", p)


if __name__ == "__main__":
    main()
