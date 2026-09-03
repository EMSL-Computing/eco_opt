"""Open-loop tape vs shrinking closed-loop on the plant–soil OCP.

Same process: closed loop should reprint the open-loop controls.
Unforecast mineral-N pulse: closed loop re-plans allocation; the tape does not.
Coarse mesh (N=8) so the demo finishes in reasonable time.
"""
from __future__ import annotations

import os
import sys

import casadi as cs
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import eco_opt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import plant_soil_nitrogen_allocation as ro

FIG_DIR = os.path.join(HERE, "figs")
N = 8
T = 120.0
DT = T / N
PULSE_T0, PULSE_TF = 30.0, 50.0
PULSE_N = 0.4
IPOPT = {"print_level": 0, "max_iter": 1500, "tol": 1e-5, "sb": "yes"}


def x0_vec():
    s = ro.STATE_INIT
    return [
        s["leafC"], s["rootC"], s["exuC"], s["SOMC"], s["micC"],
        s["rootN"], s["SOMN"], s["micN"], s["inorgN"],
    ]


def pulse_model(t, x, u):
    dx, y = ro.plant_soil_microbe_nested(t, x, u, ro.PARAMS)
    I = PULSE_N * cs.if_else(cs.logic_and(t >= PULSE_T0, t < PULSE_TF), 1.0, 0.0)
    parts = [dx[i] for i in range(8)] + [dx[8] + I]
    return cs.vertcat(*parts), y


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    ocp, sys = ro.build_ocp(terminal_time=T)
    init = {sys.x[i]: x0_vec()[i] for i in range(9)}
    init[sys.u[0]] = 1.0
    init[sys.u[1]] = 0.25
    sol = ocp.solve(
        control_intervals=N,
        polynomial_degree=3,
        collocation_points="legendre",
        initial_guess=init,
        tf_guess=T,
        ipopt_options=IPOPT,
    )
    assert sol.success, sol.return_status
    x0 = x0_vec()
    plant = eco_opt.System(states=9, controls=2, model=pulse_model, name="root_pulse")
    ol_pulse = eco_opt.simulate(plant, x0, sol.time, u=sol)
    cl_pulse = ocp.closed_loop(
        x0=x0,
        dt=DT,
        t_final=T,
        horizon="shrinking",
        plant=plant,
        polynomial_degree=3,
        collocation_points="legendre",
        ipopt_options=IPOPT,
        control_intervals=N,
        initial_guess=init,
    )
    print("root_opt  T=120 d  N=8  open tape vs closed loop on extra mineral N")
    print(f"open-loop NLP          J={sol.objective:.2f}  leafC(T)={sol.x[-1,0]:.2f}")
    print(f"tape on N pulse        leafC(T)={ol_pulse.x[-1,0]:.2f}  inorgN max={ol_pulse.x[:,8].max():.2f}")
    print(f"closed loop on N pulse leafC(T)={cl_pulse.x[-1,0]:.2f}  inorgN max={cl_pulse.x[:,8].max():.2f}")

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.4), sharex=True)
    axes[0].axvspan(PULSE_T0, PULSE_TF, color="C2", alpha=0.15, label="mineral-N pulse")
    axes[0].plot(sol.time, sol.x[:, 8], "k-", lw=2, label="open-loop inorg N")
    axes[0].plot(ol_pulse.time, ol_pulse.x[:, 8], "C0--", lw=1.5, label="tape on pulse")
    axes[0].plot(cl_pulse.time, cl_pulse.x[:, 8], "C3-", lw=2, label="closed loop")
    axes[0].set_ylabel("mineral N")
    axes[0].legend(fontsize=8)
    axes[1].step(sol.time, np.append(sol.u[:, 1], sol.u[-1, 1]), where="post", color="k", lw=2, label="tape exudation fraction")
    axes[1].step(
        cl_pulse.time, np.append(cl_pulse.u[:, 1], cl_pulse.u[-1, 1]),
        where="post", color="C3", lw=2, label="closed-loop exudation fraction",
    )
    axes[1].set_ylabel("exudation fraction")
    axes[1].set_xlabel("t (d)")
    axes[1].legend(fontsize=8)
    fig.suptitle("Plant–soil: open-loop tape vs closed loop after extra mineral N", fontsize=10)
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "open_vs_closed.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
