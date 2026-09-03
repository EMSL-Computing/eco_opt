"""Lenhart and Workman, Lab 12 (bioreactor nutrient; bang-bang).

Equivalent problem after eliminating contaminant z:

    maximize   int_0^T (K x - u) dt
    subject to xdot = G u x - D x^2,  x(0)=x0,  0 <= u <= M

Book (19.2) K=2,G=1: no injection (u≡0). (19.3) K=2.1: inject early.
z(t) = z0 exp(-K int x), recovered after the solve.
"""

import numpy as np

import eco_opt
from eco_opt import Namespace

DEFAULT = dict(K=2.0, G=1.0, D=1.0, x0=0.5, z0=0.1, M=1.0, T=2.0)
INJECT = dict(K=2.1, G=1.0, D=1.0, x0=0.5, z0=0.1, M=1.0, T=2.0)


def recover_z(sol, K, z0):
    t = sol.time
    x = sol.x[:, 0]
    integ = np.concatenate([[0.0], np.cumsum(0.5 * (x[1:] + x[:-1]) * np.diff(t))])
    return z0 * np.exp(-K * integ)


def _model(G, D):
    def f(t, x, u):
        return G * u[0] * x[0] - D * x[0] ** 2, Namespace(x=x[0], u=u[0])

    return f


def build_ocp(**kw):
    p = dict(DEFAULT)
    p.update({k: float(v) for k, v in kw.items() if v is not None})
    sys = eco_opt.System(
        states=1, controls=1, model=_model(p["G"], p["D"]), name="lab12"
    )
    ocp = eco_opt.Ocp()
    ocp.max([p["K"] * sys.y.x - sys.y.u])
    ocp.st(
        system=sys,
        constraints=[
            (0.0, "==", eco_opt.t_0(sys.t)),
            (p["T"], "==", eco_opt.t_f(sys.t)),
            (p["x0"], "==", eco_opt.t_0(sys.y.x)),
            (0.0, "<=", sys.y.u, "<=", p["M"]),
            (0.0, "<=", sys.y.x, "<=", 1e6),
        ],
    )
    return ocp, sys, p


def solve(control_intervals=60, polynomial_degree=3, **kw):
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
    print("u mean", float(sol.u.mean()))
