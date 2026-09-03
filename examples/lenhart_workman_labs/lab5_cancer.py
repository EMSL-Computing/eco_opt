"""Lenhart and Workman, Lab 5 (chemotherapy, Skipper log-kill).

    minimize   int_0^T (a N^2 + u^2) dt
    subject to Ndot = r N ln(1/N) - u delta N,  N(0)=N0 in (0,1)
               u >= 0
"""

import casadi as cs
import numpy as np

import eco_opt
from eco_opt import Namespace

DEFAULT = dict(r=0.3, a=3.0, delta=0.45, N0=0.975, T=20.0)


def untreated_N(t, r, N0):
    """Gompertz without drug: N(t) = N0 ** exp(-r t)  (carrying capacity 1)."""
    t = np.asarray(t, dtype=float)
    return N0 ** np.exp(-r * t)


def _model(r, delta):
    def f(t, x, u):
        N, ctrl = x[0], u[0]
        dx = r * N * cs.log(1.0 / N) - ctrl * delta * N
        return dx, Namespace(N=N, u=ctrl)

    return f


def build_ocp(**kw):
    p = dict(DEFAULT)
    p.update({k: float(v) for k, v in kw.items() if v is not None})
    sys = eco_opt.System(
        states=1, controls=1, model=_model(p["r"], p["delta"]), name="lab5"
    )
    ocp = eco_opt.Ocp()
    ocp.min([p["a"] * sys.y.N**2 + sys.y.u**2])
    ocp.st(
        system=sys,
        constraints=[
            (0.0, "==", eco_opt.t_0(sys.t)),
            (p["T"], "==", eco_opt.t_f(sys.t)),
            (p["N0"], "==", eco_opt.t_0(sys.y.N)),
            (1e-8, "<=", sys.y.N, "<=", 1.0 - 1e-8),
            (0.0, "<=", sys.y.u, "<=", cs.inf),
        ],
    )
    return ocp, sys, p


def solve(control_intervals=50, polynomial_degree=3, **kw):
    ocp, sys, p = build_ocp(**kw)
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        initial_guess={sys.x[0]: p["N0"], sys.u[0]: 0.1},
        tf_guess=p["T"],
    )
    return sol, sys, p


if __name__ == "__main__":
    sol, _, p = solve()
    print(sol)
    print(p)
