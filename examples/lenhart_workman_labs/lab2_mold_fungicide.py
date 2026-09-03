"""Lenhart and Workman, Lab 2 (mold and fungicide).

Problem::

    minimize   int_0^T (A x^2 + u^2) dt
    subject to xdot = r (M - x) - u x,  x(0) = x0 > 0

Default book values (6.1): r=0.3, M=10, A=1, x0=1, T=5.
A second book case (6.2) uses A=10.

No closed-form (x,u). With u ≡ 0 the state is
x(t) = M + (x0 - M) exp(-r t).
"""

import casadi as cs
import numpy as np

import eco_opt
from eco_opt import Namespace

DEFAULT = dict(r=0.3, M=10.0, A=1.0, x0=1.0, T=5.0)
CASE_A10 = dict(r=0.3, M=10.0, A=10.0, x0=1.0, T=5.0)


def untreated_state(t, r, M, x0):
    t = np.asarray(t, dtype=float)
    return M + (x0 - M) * np.exp(-r * t)


def _model(r, M):
    def mold_model(t, x, u):
        state, ctrl = x[0], u[0]
        dx = r * (M - state) - ctrl * state
        y = Namespace(x=state, u=ctrl)
        return dx, y

    return mold_model


def build_ocp(r=None, M=None, A=None, x0=None, T=None):
    p = dict(DEFAULT)
    if r is not None:
        p["r"] = float(r)
    if M is not None:
        p["M"] = float(M)
    if A is not None:
        p["A"] = float(A)
    if x0 is not None:
        p["x0"] = float(x0)
    if T is not None:
        p["T"] = float(T)

    sys = eco_opt.System(
        states=1, controls=1, model=_model(p["r"], p["M"]), name="lab2"
    )
    ocp = eco_opt.Ocp()
    ocp.min([p["A"] * sys.y.x**2 + sys.y.u**2])
    ocp.st(
        system=sys,
        constraints=[
            (0.0, "==", eco_opt.t_0(sys.t)),
            (p["T"], "==", eco_opt.t_f(sys.t)),
            (p["x0"], "==", eco_opt.t_0(sys.y.x)),
            (0.0, "<=", sys.y.x, "<=", cs.inf),
            (-cs.inf, "<=", sys.y.u, "<=", cs.inf),
        ],
    )
    return ocp, sys, p


def solve(
    r=None,
    M=None,
    A=None,
    x0=None,
    T=None,
    control_intervals=50,
    polynomial_degree=3,
):
    ocp, sys, p = build_ocp(r=r, M=M, A=A, x0=x0, T=T)
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        initial_guess={sys.x[0]: p["x0"], sys.u[0]: 0.0},
        tf_guess=p["T"],
    )
    return sol, sys, p


if __name__ == "__main__":
    sol, _sys, p = solve()
    print(sol)
    print("params", p)
    print("u[0], u[-1]:", float(sol.u[0, 0]), float(sol.u[-1, 0]))
    print("x[0], x[-1]:", float(sol.x[0, 0]), float(sol.x[-1, 0]))
