"""Lenhart and Workman, Lab 6 (fish harvesting).

    maximize   int_0^T [ A k t/(t+1) x u - u^2 ] dt
    subject to xdot = -(m + u) x,  x(0)=x0
               0 <= u <= M
"""

import casadi as cs
import numpy as np

import eco_opt
from eco_opt import Namespace

DEFAULT = dict(A=5.0, k=10.0, m=0.2, x0=0.4, M=1.0, T=10.0)


def unharvested(t, m, x0):
    return x0 * np.exp(-m * np.asarray(t, dtype=float))


def _model(m):
    def f(t, x, u):
        return -(m + u[0]) * x[0], Namespace(x=x[0], u=u[0])

    return f


def build_ocp(**kw):
    p = dict(DEFAULT)
    p.update({k: float(v) for k, v in kw.items() if v is not None})
    sys = eco_opt.System(states=1, controls=1, model=_model(p["m"]), name="lab6")
    t, x, u = sys.t, sys.y.x, sys.y.u
    mass = p["k"] * t / (t + 1.0)
    ocp = eco_opt.Ocp()
    ocp.max([p["A"] * mass * x * u - u**2])
    ocp.st(
        system=sys,
        constraints=[
            (0.0, "==", eco_opt.t_0(sys.t)),
            (p["T"], "==", eco_opt.t_f(sys.t)),
            (p["x0"], "==", eco_opt.t_0(x)),
            (0.0, "<=", u, "<=", p["M"]),
        ],
    )
    return ocp, sys, p


def solve(control_intervals=50, polynomial_degree=3, **kw):
    ocp, sys, p = build_ocp(**kw)
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        initial_guess={sys.x[0]: p["x0"], sys.u[0]: 0.0},
        tf_guess=p["T"],
    )
    return sol, sys, p


if __name__ == "__main__":
    sol, _, p = solve()
    print(sol)
    print("u[0], u max, u[-1]", float(sol.u[0, 0]), float(sol.u.max()), float(sol.u[-1, 0]))
