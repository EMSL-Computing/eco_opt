"""Open-loop tape vs shrinking closed-loop on the lignin OCP.

Same process model: closed loop reprints the open-loop v_O.
Extra litter C mid-horizon: the tape ignores it; closed loop re-plans.
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
import lignin_decomposition as ld

FIG_DIR = os.path.join(HERE, "figs")

N = 12
DT = ld.TT_DEFAULT / N
PULSE_T0, PULSE_TF = 600.0, 900.0
PULSE_I = 0.02
IPOPT = {"print_level": 0, "max_iter": 2000, "sb": "yes"}


def pulse_model(t, x, u, p):
    dx, y = ld.lignin_model(t, x, u, p)
    I = PULSE_I * cs.if_else(cs.logic_and(t >= PULSE_T0, t < PULSE_TF), 1.0, 0.0)
    return [dx[0] + I, dx[1]], y


def growth_rate(CT, Co, vo, p):
    vh_max, _mo, ro = np.asarray(p, dtype=float).ravel()
    L = Co / CT
    Dh = vh_max * np.exp(-((L / ld.A) ** ld.B)) * (CT - Co)
    Do = vo * Co
    e = ld.EMAX - ro * vo
    return e * (Dh + Do)


def cumulative_G(time, x, u, p):
    """Realized ∫G dt on the applied piecewise-constant \(v_O\)."""
    t = np.asarray(time, dtype=float)
    vo = np.asarray(u, dtype=float).reshape(-1)
    acc = 0.0
    for k in range(vo.size):
        acc += growth_rate(x[k, 0], x[k, 1], vo[k], p) * (t[k + 1] - t[k])
    return float(acc)


def closed_loop(ocp, sys, x0, plant=None, p=None):
    return ocp.closed_loop(
        x0=x0,
        dt=DT,
        t_final=ld.TT_DEFAULT,
        horizon="shrinking",
        plant=plant,
        p=p,
        polynomial_degree=3,
        collocation_points="legendre",
        ipopt_options=IPOPT,
        control_intervals=N,
    )


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    ocp, sys = ld.build_ocp(ld.VH_LOOP, ld.MO_DEFAULT, ld.RO_LOOP)
    sol = ocp.solve(
        control_intervals=N,
        polynomial_degree=3,
        collocation_points="legendre",
        tf_guess=ld.TT_DEFAULT,
        ipopt_options=IPOPT,
    )
    assert sol.success, sol.return_status
    x0 = sol.x[0].tolist()
    p = sol.p
    replay = ocp.simulate(sol)
    plant = eco_opt.System(
        states=2, controls=1, parameters=3, model=pulse_model, name="lignin_pulse"
    )
    ol_pulse = eco_opt.simulate(plant, x0, sol.time, u=sol, p=p)
    cl_same = closed_loop(ocp, sys, x0, plant=None, p=p)
    cl_pulse = closed_loop(ocp, sys, x0, plant=plant, p=p)

    J_ol = float(sol.objective)
    J_replay = cumulative_G(replay.time, replay.x, replay.u, p)
    J_cl = cumulative_G(cl_same.time, cl_same.x, cl_same.u, p)
    print("Lignin  T=2000 d  open vs closed")
    print(f"open-loop NLP          J=∫G = {J_ol:.4f}  CT(T)={sol.x[-1,0]:.4f}")
    print(f"replay u* on f         J=∫G = {J_replay:.4f}  CT(T)={replay.x[-1,0]:.4f}")
    print(f"closed loop, same f    J=∫G = {J_cl:.4f}  CT(T)={cl_same.x[-1,0]:.4f}")
    print(f"tape on extra-C pulse  CT(T)={ol_pulse.x[-1,0]:.4f}")
    print(f"closed loop on pulse   CT(T)={cl_pulse.x[-1,0]:.4f}")

    t = sol.time
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.2), sharex="col")
    axes[0, 0].plot(t, sol.x[:, 0], "k-", lw=2, label="open-loop $C_T$")
    axes[0, 0].plot(cl_same.time, cl_same.x[:, 0], "C1:", lw=2, label="closed loop")
    axes[0, 0].set_title("same process as the OCP")
    axes[0, 0].set_ylabel("$C_T$")
    axes[0, 0].legend(fontsize=8)
    axes[1, 0].step(t, np.append(sol.u[:, 0], sol.u[-1, 0]), where="post", color="k", lw=2)
    axes[1, 0].step(
        cl_same.time, np.append(cl_same.u[:, 0], cl_same.u[-1, 0]),
        where="post", color="C1", linestyle=":",
    )
    axes[1, 0].set_ylabel("$v_O$")
    axes[1, 0].set_xlabel("t (d)")
    axes[0, 1].axvspan(PULSE_T0, PULSE_TF, color="C2", alpha=0.15, label="extra C")
    axes[0, 1].plot(ol_pulse.time, ol_pulse.x[:, 0], "k-", lw=2, label="open-loop tape")
    axes[0, 1].plot(cl_pulse.time, cl_pulse.x[:, 0], "C3-", lw=2, label="closed loop")
    axes[0, 1].set_title("extra litter C in [600, 900] d")
    axes[0, 1].set_ylabel("$C_T$")
    axes[0, 1].legend(fontsize=8)
    axes[1, 1].axvspan(PULSE_T0, PULSE_TF, color="C2", alpha=0.15)
    axes[1, 1].step(
        ol_pulse.time, np.append(ol_pulse.u[:, 0], ol_pulse.u[-1, 0]),
        where="post", color="k", lw=2, label="tape",
    )
    axes[1, 1].step(
        cl_pulse.time, np.append(cl_pulse.u[:, 0], cl_pulse.u[-1, 0]),
        where="post", color="C3", lw=2, label="closed loop",
    )
    axes[1, 1].set_ylabel("$v_O$")
    axes[1, 1].set_xlabel("t (d)")
    fig.suptitle("Ligninolytic rate: open-loop tape vs shrinking closed loop", fontsize=10)
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "open_vs_closed.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
