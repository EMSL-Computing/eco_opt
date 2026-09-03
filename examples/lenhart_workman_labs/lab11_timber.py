"""Lenhart and Workman, Lab 11 (timber reinvestment; bang-bang).

    maximize   int_0^T exp(-r t) x (1-u) dt
    subject to xdot = k x u,  x(0)=x0,  0 <= u <= 1

Analytic (book): no singular arc. For r=0 the terminal no-investment
interval has length 1/k (independent of x0). With r=0, k=1, T=5 the
switch is at t*=4: u=1 on [0,4), u=0 on (4,5].
"""

import casadi as cs

import eco_opt
from eco_opt import Namespace

DEFAULT = dict(x0=100.0, k=1.0, r=0.0, T=5.0)


def analytic_switch(k, r, T):
    """Switch time for r=0 is T - 1/k, clipped to [0, T]."""
    if abs(r) < 1e-12:
        return min(max(T - 1.0 / k, 0.0), T)
    return None


def _model(k):
    def f(t, x, u):
        return k * x[0] * u[0], Namespace(x=x[0], u=u[0])

    return f


def build_ocp(**kw):
    p = dict(DEFAULT)
    p.update({k: float(v) for k, v in kw.items() if v is not None})
    sys = eco_opt.System(states=1, controls=1, model=_model(p["k"]), name="lab11")
    ocp = eco_opt.Ocp()
    ocp.max([cs.exp(-p["r"] * sys.t) * sys.y.x * (1.0 - sys.y.u)])
    ocp.st(
        system=sys,
        constraints=[
            (0.0, "==", eco_opt.t_0(sys.t)),
            (p["T"], "==", eco_opt.t_f(sys.t)),
            (p["x0"], "==", eco_opt.t_0(sys.y.x)),
            (0.0, "<=", sys.y.u, "<=", 1.0),
            (0.0, "<=", sys.y.x, "<=", 1e12),
        ],
    )
    return ocp, sys, p


def solve(control_intervals=80, polynomial_degree=3, **kw):
    ocp, sys, p = build_ocp(**kw)
    sol = ocp.solve(
        control_intervals=control_intervals,
        polynomial_degree=polynomial_degree,
        collocation_points="legendre",
        initial_guess={sys.x[0]: p["x0"], sys.u[0]: 1.0},
        tf_guess=p["T"],
    )
    return sol, sys, p


if __name__ == "__main__":
    sol, _, p = solve()
    print(sol)
    print("analytic switch", analytic_switch(p["k"], p["r"], p["T"]))
