"""Goddard rocket, maximum ascent, free final time.

Run directly:  python examples/goddard_benchmark/goddard_benchmark.py
"""

import numpy as np
import casadi as cs

import eco_opt
from eco_opt import Namespace


def goddard_model(t, x, u):
    v, h, m = x[0], x[1], x[2]
    F = u[0]

    D0 = 0.01227
    beta = 0.145e-3
    c = 2060.0
    g0 = 9.81
    r0 = 6.371e6

    D = D0 * cs.exp(-beta * h)
    F_D = cs.sign(v) * D * v ** 2
    g = g0 * (r0 / (r0 + h)) ** 2

    dv = (F * c - F_D) / m - g
    dh = v
    dm = -F
    dx = [dv, dh, dm]

    y = Namespace()
    y.rocket = Namespace(velocity=v, height=h, mass=m, fuelMassFlow=F)
    y.drag = Namespace(coefficient=D, force=F_D)
    y.gravity = g
    return dx, y


def build_ocp(control_intervals=59):
    sys = eco_opt.System(states=3, controls=1, model=goddard_model, name="goddard")
    time = sys.t
    rocket = sys.y.rocket

    m0 = 214.839
    mf = 67.9833
    Fm = 9.525515

    ocp = eco_opt.Ocp()
    ocp.max(eco_opt.t_f(rocket.height))
    ocp.st(
        system=sys,
        constraints=[
            (0, "==", eco_opt.t_0(time)),
            (0, "==", eco_opt.t_0(rocket.velocity)),
            (0, "==", eco_opt.t_0(rocket.height)),
            (m0, "==", eco_opt.t_0(rocket.mass)),
            (0, "<=", eco_opt.t_f(time), "<=", float("inf")),
            (0, "<=", rocket.velocity, "<=", float("inf")),
            (0, "<=", rocket.height, "<=", float("inf")),
            (mf, "<=", rocket.mass, "<=", m0),
            (0, "<=", rocket.fuelMassFlow, "<=", Fm),
        ],
    )
    return ocp, sys


def solve(control_intervals=59, polynomial_degree=5, tf_guess=200.0):
    ocp, sys = build_ocp(control_intervals)
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        tf_guess=tf_guess,
    )
    return sol, sys


if __name__ == "__main__":
    sol, sys = solve()
    print(sol)
    print("objective (final height):", sol.objective)
    print("final time:", sol.tf)
