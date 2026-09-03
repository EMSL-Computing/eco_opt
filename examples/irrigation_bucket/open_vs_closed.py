"""Dry open-loop irrigation vs shrinking closed-loop, with and without rain."""
from __future__ import annotations

import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import eco_opt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import irrigation_bucket as ib

FIG_DIR = os.path.join(HERE, "figs")
DT = ib.T_HORIZON / ib.N_INTERVALS
IPOPT = {"max_iter": 3000, "print_level": 0, "sb": "yes"}


def water_use(t, u):
    return float(np.sum(np.asarray(u).ravel() * np.diff(t)))


def closed_loop(plant=None):
    ocp, sys = ib.build_ocp()
    return ocp.closed_loop(
        x0=[ib.S0],
        dt=DT,
        t_final=ib.T_HORIZON,
        horizon="shrinking",
        plant=plant,
        polynomial_degree=ib.DEGREE,
        collocation_points="radau",
        ipopt_options=IPOPT,
        initial_guess={sys.u[0]: ib.E0},
    )


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    sol, sys, ocp = ib.solve()
    assert sol.success, sol.return_status
    replay = ocp.simulate(sol)
    plant = ib.build_rain_plant()
    ol_rain = eco_opt.simulate(plant, [ib.S0], sol.time, u=sol)
    cl_dry = closed_loop(plant=None)
    cl_rain = closed_loop(plant=plant)

    J_ol = float(sol.objective)
    J_cl = water_use(cl_dry.time, cl_dry.u[:, 0])
    J_ol_rain = water_use(ol_rain.time, ol_rain.u[:, 0])
    J_cl_rain = water_use(cl_rain.time, cl_rain.u[:, 0])
    du = np.max(np.abs(cl_dry.u[:, 0] - sol.u[:, 0]))

    print("Irrigation bucket  T=14 d  dry OCP (P=0, sinusoidal E)")
    print(f"open-loop NLP              ∫u = {J_ol:.3f}   S(T)={sol.x[-1,0]:.3f}  min S={sol.x[:,0].min():.3f}")
    print(f"replay u* on dry plant     ∫u = {water_use(replay.time, replay.u[:,0]):.3f}   S(T)={replay.x[-1,0]:.3f}")
    print(f"closed loop, dry plant     ∫u = {J_cl:.3f}   S(T)={cl_dry.x[-1,0]:.3f}  max|u-u*|={du:.4f}")
    print(f"open-loop tape on rain     ∫u = {J_ol_rain:.3f}   S(T)={ol_rain.x[-1,0]:.3f}  max S={ol_rain.x[:,0].max():.3f}")
    print(f"closed loop on rain        ∫u = {J_cl_rain:.3f}   S(T)={cl_rain.x[-1,0]:.3f}  max S={cl_rain.x[:,0].max():.3f}")

    t_ol = sol.time
    u_ol = np.append(sol.u[:, 0], sol.u[-1, 0])
    E = ib.E0 * (1.0 + ib.E_AMP * np.sin(2.0 * np.pi * t_ol))

    fig, axes = plt.subplots(3, 2, figsize=(8.4, 7.4), sharex="col")

    ax = axes[0, 0]
    ax.plot(t_ol, sol.x[:, 0], "k-", lw=2, label="open-loop $S^*$")
    ax.plot(replay.time, replay.x[:, 0], "C0--", lw=1.4, label="replay $u^*$")
    ax.plot(cl_dry.time, cl_dry.x[:, 0], "C1:", lw=2, label="closed loop")
    ax.axhline(ib.S_MIN, color="0.5", ls="--", lw=0.8)
    ax.set_ylabel("S (mm)")
    ax.set_title("dry plant = OCP")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1, 0]
    ax.step(t_ol, u_ol, where="post", color="k", lw=2, label="open-loop $u^*$")
    ax.step(
        cl_dry.time,
        np.append(cl_dry.u[:, 0], cl_dry.u[-1, 0]),
        where="post",
        color="C1",
        lw=1.5,
        linestyle=":",
        label="closed loop",
    )
    ax.set_ylabel("u (mm d$^{-1}$)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[2, 0]
    ax.plot(t_ol, E, color="C2", lw=1.4)
    ax.set_ylabel("E(t) (mm d$^{-1}$)")
    ax.set_xlabel("t (d)")
    ax.grid(alpha=0.3)

    ax = axes[0, 1]
    ax.axvspan(ib.PULSE_T0, ib.PULSE_TF, color="C0", alpha=0.15, label="rain")
    ax.plot(ol_rain.time, ol_rain.x[:, 0], "k-", lw=2, label="open-loop tape")
    ax.plot(cl_rain.time, cl_rain.x[:, 0], "C3-", lw=2, label="closed loop")
    ax.axhline(ib.S_MIN, color="0.5", ls="--", lw=0.8)
    ax.axhline(ib.S_MAX, color="0.5", ls=":", lw=0.8)
    ax.set_ylabel("S (mm)")
    ax.set_title("plant rains on days 5–7")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    ax.axvspan(ib.PULSE_T0, ib.PULSE_TF, color="C0", alpha=0.15)
    ax.step(
        ol_rain.time,
        np.append(ol_rain.u[:, 0], ol_rain.u[-1, 0]),
        where="post",
        color="k",
        lw=2,
        label="open-loop tape",
    )
    ax.step(
        cl_rain.time,
        np.append(cl_rain.u[:, 0], cl_rain.u[-1, 0]),
        where="post",
        color="C3",
        lw=2,
        label="closed loop",
    )
    ax.set_ylabel("u (mm d$^{-1}$)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[2, 1]
    ax.axvspan(ib.PULSE_T0, ib.PULSE_TF, color="C0", alpha=0.15)
    t_fine = np.linspace(0, ib.T_HORIZON, 400)
    P = np.where((t_fine >= ib.PULSE_T0) & (t_fine < ib.PULSE_TF), ib.P_RAIN, 0.0)
    ax.plot(t_fine, P, color="C0", lw=1.4)
    ax.set_ylabel("P(t) (mm d$^{-1}$)")
    ax.set_xlabel("t (d)")
    ax.grid(alpha=0.3)

    fig.suptitle(
        "Irrigation bucket: open loop vs shrinking closed loop\n"
        f"dry: $\\int u$={J_ol:.1f} (ol) / {J_cl:.1f} (cl)  |  "
        f"rain plant: tape {J_ol_rain:.1f}, closed {J_cl_rain:.1f} mm",
        fontsize=10,
    )
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "open_vs_closed.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
