"""Open-loop tape vs shrinking closed-loop on the microbial growth OCP.

Same plant: closed loop reprints the open-loop uptake.
Plant with an unexpected C pulse: closed loop sees the extra substrate
and re-plans; the open-loop tape does not.
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
from eco_opt import Namespace

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import microbial_growth as mg

FIG_DIR = os.path.join(HERE, "figs")

T = 1.50
N = 15
DT = T / N
PULSE_T0, PULSE_TF = 0.40, 0.70
PULSE_I = 1.2
IPOPT = {"print_level": 0, "max_iter": 3000, "sb": "yes"}


def growth(u):
    return mg.growth_rate(u)


def ocp_model(t, x, u):
    C, uptake = x[0], u[0]
    G = growth(uptake)
    dC = -uptake - mg.GAMMA * C + mg.MU * G
    return dC, Namespace(C=C, u=uptake, G=G)


def pulse_plant(t, x, u):
    C, uptake = x[0], u[0]
    G = growth(uptake)
    I = PULSE_I * cs.if_else(cs.logic_and(t >= PULSE_T0, t < PULSE_TF), 1.0, 0.0)
    dC = -uptake - mg.GAMMA * C + mg.MU * G + I
    return dC, Namespace(C=C, u=uptake, G=G, I=I)


def build_fixed_T():
    sys = eco_opt.System(states=1, controls=1, model=ocp_model, name="microbe")
    ocp = eco_opt.Ocp()
    ocp.max([sys.y.G])
    ocp.st(
        system=sys,
        constraints=[
            (0.0, "<=", sys.x[0], "<=", 2.0),
            (mg.U_MIN, "<=", sys.u[0], "<=", mg.U_MAX),
            (0.0, "==", eco_opt.t_0(sys.t)),
            (T, "==", eco_opt.t_f(sys.t)),
            (mg.C0, "==", eco_opt.t_0(sys.x[0])),
        ],
    )
    return ocp, sys


def cumulative_G(t, u):
    """∫ G(u) dt with piecewise-constant u on t[:-1]."""
    dt = np.diff(t)
    g = growth(np.asarray(u).ravel())
    return float(np.sum(g * dt))


def run_open_loop():
    ocp, sys = build_fixed_T()
    sol = ocp.solve(
        control_intervals=N,
        polynomial_degree=3,
        collocation_points="radau",
        initial_guess={sys.x[0]: mg.C0, sys.u[0]: 0.5},
        tf_guess=T,
        ipopt_options=IPOPT,
    )
    return ocp, sys, sol


def run_closed(plant=None):
    ocp, sys = build_fixed_T()
    return ocp.closed_loop(
        x0=[mg.C0],
        dt=DT,
        t_final=T,
        horizon="shrinking",
        plant=plant,
        polynomial_degree=3,
        collocation_points="radau",
        ipopt_options=IPOPT,
        initial_guess={sys.u[0]: 0.5},
    )


def print_table(title, t, C, u, n=8):
    print(f"\n{title}")
    print(f"{'t':>6} {'C':>9} {'u':>9} {'G(u)':>9}")
    idx = np.unique(np.round(np.linspace(0, len(u) - 1, n)).astype(int))
    for k in idx:
        print(f" {t[k]:5.2f} {C[k]:9.4f} {u[k]:9.4f} {float(growth(u[k])):9.4f}")
    print(f" {t[-1]:5.2f} {C[-1]:9.4f} {'(end)':>9}")


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    ocp, sys, sol = run_open_loop()
    assert sol.success, sol.return_status
    replay = ocp.simulate(sol)
    plant = eco_opt.System(states=1, controls=1, model=pulse_plant, name="pulse")
    ol_on_pulse = eco_opt.simulate(plant, [mg.C0], sol.time, u=sol)
    cl_same = run_closed(plant=None)
    cl_pulse = run_closed(plant=plant)

    J_ol = float(sol.objective)
    J_replay = cumulative_G(replay.time, replay.u[:, 0])
    J_cl = cumulative_G(cl_same.time, cl_same.u[:, 0])
    J_ol_pulse = cumulative_G(ol_on_pulse.time, ol_on_pulse.u[:, 0])
    J_cl_pulse = cumulative_G(cl_pulse.time, cl_pulse.u[:, 0])

    print("Microbial growth, fixed T=1.5 (Manzoni-style C, u, G)")
    print(f"open-loop NLP     J=∫G = {J_ol:.4f}   C(T)={sol.x[-1,0]:.4f}  status={sol.return_status}")
    print(f"replay u* on f    J=∫G = {J_replay:.4f}   C(T)={replay.x[-1,0]:.4f}  (same model)")
    print(
        f"closed loop, same f J=∫G = {J_cl:.4f}   C(T)={cl_same.x[-1,0]:.4f}  "
        f"max|u-u*|={np.max(np.abs(cl_same.u[:,0]-sol.u[:,0])):.4f}"
    )
    print(
        f"open-loop tape on pulsed plant  J=∫G = {J_ol_pulse:.4f}   C(T)={ol_on_pulse.x[-1,0]:.4f}"
    )
    print(
        f"closed loop on pulsed plant     J=∫G = {J_cl_pulse:.4f}   C(T)={cl_pulse.x[-1,0]:.4f}"
    )

    print_table("OPEN-LOOP u*(t)", sol.time, sol.x[:, 0], sol.u[:, 0])
    print_table("CLOSED-LOOP, same plant", cl_same.time, cl_same.x[:, 0], cl_same.u[:, 0])
    print_table("OPEN-LOOP tape on pulse plant", ol_on_pulse.time, ol_on_pulse.x[:, 0], ol_on_pulse.u[:, 0])
    print_table("CLOSED-LOOP on pulse plant", cl_pulse.time, cl_pulse.x[:, 0], cl_pulse.u[:, 0])

    t_ol = sol.time
    u_ol = np.append(sol.u[:, 0], sol.u[-1, 0])
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.2), sharex="col")

    ax = axes[0, 0]
    ax.plot(t_ol, sol.x[:, 0], "k-", lw=2, label="open-loop $C^*$")
    ax.plot(replay.time, replay.x[:, 0], "C0--", lw=1.5, label="replay $u^*$")
    ax.plot(cl_same.time, cl_same.x[:, 0], "C1:", lw=2, label="closed loop")
    ax.set_ylabel("C (substrate)")
    ax.set_title("same plant as the OCP")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1, 0]
    ax.step(t_ol, u_ol, where="post", color="k", lw=2, label="open-loop $u^*$")
    ax.step(cl_same.time, np.append(cl_same.u[:, 0], cl_same.u[-1, 0]), where="post",
            color="C1", lw=1.5, linestyle=":", label="closed loop")
    ax.set_ylabel("u (uptake)")
    ax.set_xlabel("t")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[0, 1]
    ax.axvspan(PULSE_T0, PULSE_TF, color="C2", alpha=0.15, label="C pulse")
    ax.plot(ol_on_pulse.time, ol_on_pulse.x[:, 0], "k-", lw=2, label="open-loop tape")
    ax.plot(cl_pulse.time, cl_pulse.x[:, 0], "C3-", lw=2, label="closed loop")
    ax.set_ylabel("C (substrate)")
    ax.set_title("plant gets extra C in [0.4, 0.7]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    ax.axvspan(PULSE_T0, PULSE_TF, color="C2", alpha=0.15)
    ax.step(ol_on_pulse.time, np.append(ol_on_pulse.u[:, 0], ol_on_pulse.u[-1, 0]),
            where="post", color="k", lw=2, label="open-loop tape")
    ax.step(cl_pulse.time, np.append(cl_pulse.u[:, 0], cl_pulse.u[-1, 0]),
            where="post", color="C3", lw=2, label="closed loop")
    ax.set_ylabel("u (uptake)")
    ax.set_xlabel("t")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.suptitle(
        "Microbial uptake: open loop vs shrinking closed loop\n"
        f"same $f$: $J_\\mathrm{{ol}}$={J_ol:.3f}, $J_\\mathrm{{cl}}$={J_cl:.3f}  |  "
        f"pulse plant: tape $J$={J_ol_pulse:.3f}, closed $J$={J_cl_pulse:.3f}",
        fontsize=10,
    )
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "open_vs_closed.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
