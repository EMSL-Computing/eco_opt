"""Lenhart and Workman, Lab 4 (Lab 1 with bounds on u).

    maximize   int_0^1 (A x - u^2) dt
    subject to xdot = -0.5 x^2 + C u,  x(0)=x0 > -2
               M1 <= u <= M2,  A >= 0
"""

import casadi as cs

import eco_opt
from eco_opt import Namespace

T_HORIZON = 1.0
# Book (9.1): bounds contain the Lab 1 unconstrained solution.
DEFAULT = dict(A=1.0, C=4.0, x0=1.0, M1=-1.0, M2=2.0)


def _model(C):
    def f(t, x, u):
        y = Namespace(x=x[0], u=u[0])
        return -0.5 * x[0] ** 2 + C * u[0], y

    return f


def build_ocp(A=None, C=None, x0=None, M1=None, M2=None):
    p = dict(DEFAULT)
    for k, v in dict(A=A, C=C, x0=x0, M1=M1, M2=M2).items():
        if v is not None:
            p[k] = float(v)
    sys = eco_opt.System(states=1, controls=1, model=_model(p["C"]), name="lab4")
    ocp = eco_opt.Ocp()
    ocp.max([p["A"] * sys.y.x - sys.y.u**2])
    ocp.st(
        system=sys,
        constraints=[
            (0.0, "==", eco_opt.t_0(sys.t)),
            (T_HORIZON, "==", eco_opt.t_f(sys.t)),
            (p["x0"], "==", eco_opt.t_0(sys.y.x)),
            (p["M1"], "<=", sys.y.u, "<=", p["M2"]),
        ],
    )
    return ocp, sys, p


def solve(control_intervals=40, polynomial_degree=3, **kw):
    ocp, sys, p = build_ocp(**kw)
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        initial_guess={sys.x[0]: p["x0"], sys.u[0]: 0.5 * (p["M1"] + p["M2"])},
        tf_guess=T_HORIZON,
    )
    return sol, sys, p


if __name__ == "__main__":
    sol, _, p = solve()
    print(sol)
    print(p)
